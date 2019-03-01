import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TIntArrayList;
import htsjdk.samtools.*;
import org.jdom2.Document;
import org.jdom2.Element;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

public class DebugPrinter extends Constants{

  Aligner aligner;
  KMerCounter counter;
  ReferenceAlignment refAlignment;
  private Document document;
  private String pathToLongestPaths;
  private String pathToMergedContigs;
  private String pathToMergedContigsBam;
  private String pathToContigBlasts;
  private String pathToContigsCov;
  private String pathToSelectedReferences;
  private String pathToClusters;
  private String pathToCutPaths;
  private String pathToConsensuses;
  private int minContigLen;
  private int minReadNum;
  String outPath;

  private static DebugPrinter instance;

  private DebugPrinter(Document document){
    this.document = document;
    Element referenceSelectorElement = document.getRootElement().getChild("ReferenceSelector");
    if(Boolean.parseBoolean(referenceSelectorElement.getChildText("Enabled"))){
      refAlignment = ReferenceAlignment.getInstance(document);
    }
    outPath = document.getRootElement().getChildText("OutPath");
    pathToLongestPaths = outPath + "longestPaths.txt";
    pathToMergedContigs = outPath + "mergedContigs.txt";
    pathToMergedContigsBam = outPath + "mergedContigs.bam";
    pathToContigBlasts = outPath + "references.txt";
    pathToContigsCov = outPath + "mergedContigs_coverages.txt";
    pathToSelectedReferences = outPath + "selected_references.txt";
    pathToClusters = outPath + "clusters_reads.txt";
    pathToCutPaths = outPath + "cut_contigs.txt";
    pathToConsensuses = outPath + "consensuses.fasta";
    minContigLen = Integer.parseInt(referenceSelectorElement.getChildText("MinContigLength"));
    minReadNum = Integer.parseInt(
        referenceSelectorElement.getChild("Graph").getChildText("MinReadNumber"));
    aligner = new SmithWatermanGotoh(document);
    counter = KMerCounter.getInstance(document);
    try{
      File outFolder = new File(outPath);
      if(!outFolder.exists()){
        outFolder.mkdirs();
      }
    }catch(Exception e){
      e.printStackTrace();
    }
  }

  static DebugPrinter getInstance(Document document){
    String outPath = document.getRootElement().getChildText("OutPath");
    if(instance == null){
      instance = new DebugPrinter(document);
    }else if(!instance.outPath.equals(outPath)){
      instance = new DebugPrinter(document);
    }
    return instance;
  }

  private static void setCigar(MappedRead read, SAMRecord record){
    Cigar cigar = new Cigar();
    if(read.aln.start1 > 0){
      CigarElement se = new CigarElement(read.aln.start1, CigarOperator.S);
      cigar.add(se);
    }
    CigarOperator operator = CigarOperator.M;
    int len = 0;
    int mNum = 0;
    for(int i = 0; i < read.aln.length; i++){
      byte c = read.aln.sequence1[i];
      if(c == GAP){
        if(operator == CigarOperator.D){
          len ++;
        }else{
          CigarElement element = new CigarElement(len, operator);
          cigar.add(element);
          operator = CigarOperator.D;
          len = 1;
        }
      }else if(c < 'a'){
        if(operator == CigarOperator.M){
          len ++;
        }else{
          CigarElement element = new CigarElement(len, operator);
          cigar.add(element);
          operator = CigarOperator.M;
          len = 1;
        }
        mNum ++;
      }else{
        if(operator == CigarOperator.I){
          len ++;
        }else{
          CigarElement element = new CigarElement(len, operator);
          cigar.add(element);
          operator = CigarOperator.I;
          len = 1;
        }
      }
    }
    CigarElement element = new CigarElement(len, operator);
    cigar.add(element);
    if(read.aln.end1 < read.seq.length){
      CigarElement se =
          new CigarElement(read.seq.length - read.aln.end1, CigarOperator.S);
      cigar.add(se);
    }
    record.setAttribute(ReservedTagConstants.XN, mNum - read.aln.identity);
    record.setCigar(cigar);
  }

  private static void setRecordFieldsMapped(MappedRead read, SAMRecord record){
    int start1 = read.start + read.aln.start2 + 1;
    record.setAlignmentStart(start1);
    record.setReadName(read.name);
    record.setReadString(new String(read.seq));
    setCigar(read, record);
    record.setMappingQuality(255);
    record.setReferenceIndex(0);
    record.setReadNegativeStrandFlag(read.reverse == 1);
    record.setReadPairedFlag(false);
    record.setReadUnmappedFlag(false);
  }

  private class ContigMapper implements Runnable{

    public MappedRead mappedContig;
    private Cluster cluster;
    private Reference ref;

    public ContigMapper(Cluster cluster, Reference ref){
      this.cluster = cluster;
      this.ref = ref;
    }

    public HashMap<String, int[]> buildIndex(String seq, int K){
      HashMap<String, TIntArrayList> indexTemp = new HashMap<>();
      int centroidLen = seq.length();
      for(int j = 0; j <= centroidLen - K; j++){
        String s = seq.substring(j, j + K);
        TIntArrayList l = indexTemp.get(s);
        if(l == null){
          l = new TIntArrayList();
          indexTemp.put(s, l);
        }
        l.add(j);
      }
      HashMap<String, int[]> centroidIndex = new HashMap<>();
      for(Map.Entry<String, TIntArrayList> entry: indexTemp.entrySet()){
        centroidIndex.put(entry.getKey(), entry.getValue().toArray());
        entry.setValue(null);
      }
      return centroidIndex;
    }

    @Override
    public void run(){
      TByteArrayList builder = new TByteArrayList();
      for(int j = 0; j < cluster.contig.seq.length; j++){
        byte c = cluster.contig.seq[j];
        if(c != GAP){
          builder.add(c);
        }
      }
      HashMap<String, int[]> index = buildIndex(ref.seq, counter.K);
      byte[] seq = builder.toArray();
      mappedContig = counter.mapReadToRegion(seq, index, 0, ref.seq.length());
      mappedContig.seq = seq;
      mappedContig.name = cluster.name;
      mappedContig.aln = aligner.align(mappedContig.seq,
          ref.seq.substring(mappedContig.start, mappedContig.end).getBytes());
    }
  }

  private void printContigsToBAM(Collection<Cluster> vertices, String refID, String path) throws InterruptedException, IOException{

    int threadNum = Runtime.getRuntime().availableProcessors();
    ArrayList<ContigMapper> contigMappingTasks = new ArrayList<>();
    ExecutorService executor = Executors.newFixedThreadPool(threadNum);
    Reference refSeq = refAlignment.refSeqs.get(refID);
    int genomeLen = refSeq.seq.length();
    for(Cluster cluster : vertices){
      if(cluster.readNum > minReadNum){
        ContigMapper task = new ContigMapper(cluster, refSeq);
        contigMappingTasks.add(task);
        executor.execute(task);
      }
    }
    executor.shutdown();
    executor.awaitTermination(10, TimeUnit.HOURS);

    SAMFileHeader header = new SAMFileHeader();
    SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(header, true, new File(path));
    header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
    SAMSequenceRecord samSequenceRecord = new SAMSequenceRecord(refSeq.name, genomeLen);
    header.addSequence(samSequenceRecord);
    ArrayList<SAMRecord> records = new ArrayList<>();
    for(ContigMapper task: contigMappingTasks){
      SAMRecord record = new SAMRecord(header);
      setRecordFieldsMapped(task.mappedContig, record);
      records.add(record);
    }
    Collections.sort(records, new SamRecordComparator());
    for(SAMRecord record: records){
      writer.addAlignment(record);
    }
    writer.close();
    Process p = Runtime.getRuntime().exec("samtools index " + path);
    p.waitFor();
  }

  void printContigsRaw(LinkedList<Cluster>[] clusters, String refID){
    try{
      for(LinkedList<Cluster> clusterList : clusters){
        if(clusterList.size() > 0){
          Cluster cluster = clusterList.getFirst();
          printContigsToBAM(clusterList, refID, outPath + "contigs_" +
              cluster.start + "_" + cluster.end + "_" + refID + ".bam");
        }
      }
    }catch(Exception e){
      e.printStackTrace();
    }
  }

  private class PathMapper implements Runnable{

    public MappedRead mappedContig;
    public Path path;

    @Override
    public void run(){
      Reference hxb2 = refAlignment.refAlns.get(0);
      mappedContig = new MappedRead(0, hxb2.seq.length(), 0);
      Alignment aln = new Alignment();
      mappedContig.aln = aln;
      mappedContig.name = "Path_" + Integer.toString(path.id);
      TByteArrayList alnSeq = new TByteArrayList();
      TByteArrayList builder = new TByteArrayList();
      int firstMatch = -1;
      int lastMatch = 0;
      for(int i = path.contigAlnStart, j = 0; i < path.contigAlnEnd; i++, j++){
        if(path.contigAln[j] == GAP){
          if(hxb2.aln.charAt(i) != GAP){
            alnSeq.add(GAP);
          }
        }else{
          builder.add(path.contigAln[j]);
          if(hxb2.aln.charAt(i) == GAP){
            alnSeq.add(lower[path.contigAln[j]]);
          }else{
            alnSeq.add(path.contigAln[j]);
            if(firstMatch == -1){
              firstMatch = i;
            }
            lastMatch = i;
          }
        }
      }
      mappedContig.seq = builder.toArray();
      int leftSoftClip = 0;
      for(int i = 0; i < alnSeq.size(); i++){
        if(alnSeq.getQuick(i) >= 'a'){
          leftSoftClip ++;
        }else{
          break;
        }
      }
      aln.start1 = leftSoftClip;
      int rightSoftClip = 0;
      for(int i = alnSeq.size() - 1; i >= 0; i--){
        if(alnSeq.getQuick(i) >= 'a'){
          rightSoftClip ++;
        }else{
          break;
        }
      }
      aln.end1 = mappedContig.seq.length - rightSoftClip;
      aln.start2 = hxb2.alnToSeqStart(firstMatch);
      aln.end2 = hxb2.alnToSeqEnd(lastMatch);
      aln.sequence1 = new byte[alnSeq.size() - leftSoftClip - rightSoftClip];
      aln.length = aln.sequence1.length;
      alnSeq.toArray(aln.sequence1, leftSoftClip, 0, aln.sequence1.length);
    }
  }

  private static class SamRecordComparator implements Comparator<SAMRecord>{

    @Override
    public int compare(SAMRecord o1, SAMRecord o2){
      return o1.getAlignmentStart() - o2.getAlignmentStart();
    }

  }

  void printLongestPaths(Collection<Path> longestPaths) throws IOException, InterruptedException{
    BufferedWriter writer = new BufferedWriter(new FileWriter(pathToLongestPaths));
    for(Path path: longestPaths){
      writer.write(">> Path " + Integer.toString(path.id) + " clusterNum = " + path.size() +
          " totalWeight = " + path.weight + "\n");
      Iterator<Cluster> iter = path.iterator();
      while(iter.hasNext()){
        Cluster vertex = iter.next();
        writer.write("> " + vertex.name + " readNum = " + vertex.readNum + "\n");
        writer.write(new String(vertex.contig.seq));
        writer.newLine();
        StringBuilder refIDs = new StringBuilder();
        for(RefSubseq refSubseq: vertex.refSubseqs){
          refIDs.append(refAlignment.refAlns.get(refSubseq.refID).name).append(' ');
        }
        refIDs.append('\n');
        writer.write(refIDs.toString());
      }
    }
    writer.close();
    writer = new BufferedWriter(new FileWriter(pathToMergedContigs));
    for(Path path: longestPaths){
      if(path.contig.length >= minContigLen){
        writer.write("> Path " + Integer.toString(path.id) + " contigLen = " + path.contig.length + "\n");
        writer.write(new String(path.contig));
        writer.newLine();
      }
    }
    writer.close();

    int threadNum = Runtime.getRuntime().availableProcessors();
    ArrayList<PathMapper> pathMappingTasks = new ArrayList<>();
    ExecutorService executor = Executors.newFixedThreadPool(threadNum);
    for(Path path: longestPaths){
      if(path.contig.length >= minContigLen){
        PathMapper task = new PathMapper();
        task.path = path;
        pathMappingTasks.add(task);
        executor.execute(task);
      }
    }
    executor.shutdown();
    executor.awaitTermination(10, TimeUnit.HOURS);

    SAMFileHeader header = new SAMFileHeader();
    SAMFileWriter writerBam = new SAMFileWriterFactory().makeBAMWriter(header, true, new File(pathToMergedContigsBam));
    header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
    ArrayList<SAMRecord> records = new ArrayList<>();
    int genomeLen = refAlignment.refAlns.get(0).seq.length();
    SAMSequenceRecord samSequenceRecord = new SAMSequenceRecord("HXB2", genomeLen);
    header.addSequence(samSequenceRecord);
    writer = new BufferedWriter(new FileWriter(pathToContigsCov));
    for(PathMapper task: pathMappingTasks){
      Path path = task.path;
      SAMRecord record = new SAMRecord(header);
      setRecordFieldsMapped(task.mappedContig, record);
      records.add(record);
      writer.write("> Path " + Integer.toString(path.id) + " contigLen = " + path.contig.length + "\n");
      Alignment aln = task.mappedContig.aln;
      int genPos = task.mappedContig.start + aln.start2;
      for(int j = aln.start1, alnPos = 0; alnPos < aln.length; alnPos++){
        if(aln.sequence1[alnPos] == GAP){
          writer.write(Integer.toString(genPos) + " " + Integer.toString(path.coverage[j]) + "\n");
          genPos++;
        }else{
          if(aln.sequence1[alnPos] < 'a'){
            writer.write(Integer.toString(genPos) + " " + Integer.toString(path.coverage[j]) + "\n");
            genPos++;
            j++;
          }else{
            j++;
          }
        }
      }
      writer.newLine();
    }
    Collections.sort(records, new SamRecordComparator());
    for(SAMRecord record: records){
      writerBam.addAlignment(record);
    }
    writerBam.close();
    writer.close();
  }

  protected void printVertices(String pathToClusters, Collection<Cluster> vertices) throws IOException{
    BufferedWriter writer = new BufferedWriter(new FileWriter(pathToClusters));
    for(Cluster vertex: vertices){
      writer.write(vertex.toString());
    }
    writer.close();
  }

  private void printSelectedRef(HashMap<String, int[]> selectedRef, ArrayList<Path> longestPaths) throws IOException, InterruptedException{
    BufferedWriter writer = new BufferedWriter(new FileWriter(pathToSelectedReferences));
    for(Map.Entry<String, int[]> entry: selectedRef.entrySet()){
      Reference seq = refAlignment.refSeqs.get(entry.getKey());
      int[] pathIDs = entry.getValue();
      writer.write(seq.name);
      for(int pathID: pathIDs){
        writer.write(" " + Integer.toString(pathID));
      }
      writer.newLine();
    }
    writer.close();
  }

  void printCutContigs(ArrayList<Path> paths, int[] selectedRefIDs){
    try{
      BufferedWriter writer = new BufferedWriter(new FileWriter(pathToCutPaths));
      for(Path path: paths){
        writer.write(Integer.toString(path.id));
        writer.newLine();
        if(path.correctContigs == null){
          writer.write("Short path\n");
          continue;
        }
        for(int i = 0; i < path.correctContigs.length; i++){
          writer.write(refAlignment.refAlns.get(selectedRefIDs[i]).name);
          writer.newLine();
          for(Contig contig : path.correctContigs[i]){
            writer.write(new String(contig.seqAln));
            writer.newLine();
          }
        }
      }
      writer.close();
    }catch(IOException e){
      e.printStackTrace();
    }
  }

  void printConsensuses(int[] selectedRefIDs){
    try{
      BufferedWriter writer = new BufferedWriter(new FileWriter(pathToConsensuses));
      for(int refID: selectedRefIDs){
        Reference seq = refAlignment.refAlns.get(refID);
        writer.write(">" + seq.name + "\n");
        writer.write(seq.consensus);
        writer.newLine();
      }
      writer.close();
    }catch(Exception e){
      e.printStackTrace();
    }
  }

  void printClustersToFiles(LinkedList<Cluster>[] clusters, ArrayList<MappedRead> reads) throws IOException{
    for(LinkedList<Cluster> clusterList : clusters){
      for(Cluster cluster : clusterList){
        BufferedWriter writer = new BufferedWriter(new FileWriter(outPath + cluster.name + ".cluster_aln"));
        writer.write(">Centroid\n");
        writer.write(cluster.centroidSeq + "\n");
        writer.write(">Consensus\n");
        if(cluster.contig == null){
          writer.write("null\n");
        }else{
          writer.write(new String(cluster.contig.seq) + "\n");
        }
        for(TIntObjectIterator<MappedRead> iter = cluster.mReads.iterator(); iter.hasNext();){
          iter.advance();
          MappedRead read = iter.value();
          writer.write(iter.key() + "\n");
          for(int i = 0; i < read.start + read.aln.start2; i++){
            writer.write('-');
          }
          writer.write(new String(read.aln.sequence1));
          for(int i = read.start + read.aln.end2; i < cluster.centroidSeq.length(); i++){
            writer.write('-');
          }
          writer.newLine();
        }
        writer.close();
        writer = new BufferedWriter(new FileWriter(outPath + cluster.name + ".cluster"));
        writer.write(">Centroid\n");
        writer.write(cluster.centroidSeq + "\n");
        writer.write(">Consensus\n");
        if(cluster.contig == null){
          writer.write("null\n");
        }else{
          writer.write(new String(cluster.contig.seq) + "\n");
        }
        for(TIntObjectIterator<MappedRead> iter = cluster.mReads.iterator(); iter.hasNext();){
          iter.advance();
          MappedRead read = iter.value();
          writer.write(iter.key() + "\n");
          writer.write(new String(read.seq) + "\n");
        }
        writer.close();
      }
    }
  }


}
