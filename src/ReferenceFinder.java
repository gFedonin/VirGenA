import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TObjectIntHashMap;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;
import org.jdom2.Document;
import org.jdom2.Element;
import org.jdom2.input.SAXBuilder;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Created by Геннадий on 03.03.2015.
 */
public class ReferenceFinder extends Constants{

  private int fileNum;
  private MappedData mappedData;
  private MapperMSA mapper;
  private MinRefSetSelector minRefSetSelector;
  Logger logger;
  private int minLen;

  private TIntArrayList[] readGroups;
  private LinkedList<Cluster>[] clusters;
  private Graph graph;
  private LongestPathFinder longestPathFinder;
  private ArrayList<Path> longestPaths;

  private String pathToUclust;
  String outPath;
  private int windowLen;
  private int winStep;
  private Document config;

  private ReferenceAlignment refAlignment;
  private ArrayList<RefSubseq>[] refSubseqs;
  KMerCounter counter;
  Aligner aligner;

  private float uclustIdentity;
  //private float identityThreshold = 0.99f;
  private int batchSize;
  private int threadNum;

  boolean debug;
  private DebugPrinter debugPrinter;

  ReferenceFinder(Document document){
    config = document;
    Element rsElement = document.getRootElement().getChild("ReferenceSelector");
    debug = Boolean.parseBoolean(rsElement.getChildText("Debug"));
    if(debug){
      debugPrinter = DebugPrinter.getInstance(document);
    }
    uclustIdentity = Float.parseFloat(rsElement.getChildText("UclustIdentity"));
    outPath = document.getRootElement().getChildText("OutPath");
    pathToUclust = rsElement.getChildText("PathToUsearch");
    minLen = Integer.parseInt(rsElement.getChildText("MinReadLength"));
    refAlignment = ReferenceAlignment.getInstance(document);
    logger = Logger.getInstance(document);
    mapper = new MapperMSA(document);
    minRefSetSelector = new MinRefSetSelector(document);
    longestPathFinder = new LongestPathFinder(document);
    counter = KMerCounter.getInstance(document);
    aligner = new SmithWatermanGotoh(document);
    threadNum = Integer.parseInt(document.getRootElement().getChildText("ThreadNumber"));
    if(threadNum == -1){
      threadNum = Runtime.getRuntime().availableProcessors();
    }
    batchSize = Integer.parseInt(document.getRootElement().getChildText("BatchSize"));
  }

  private void groupReads() throws IOException{
    readGroups = new TIntArrayList[fileNum];
    for(int i = 0; i < readGroups.length; i++){
      readGroups[i] = new TIntArrayList();
    }
    refSubseqs = new ArrayList[fileNum];
    int readNumAfterLenFilter = 0;
    for(int i = 0; i < mappedData.mappedReads.size(); i++){
      MappedRead mappedRead = mappedData.mappedReads.get(i);
      if(mappedRead.seq.length < minLen || mappedData.isConcordant[i] == 0){//
        continue;
      }
      readNumAfterLenFilter ++;
      int j = mappedRead.start/winStep;
      while(j >= 0 && mappedRead.end <= j*winStep + windowLen){
        if(j < readGroups.length){
          readGroups[j].add(i);
        }
        j--;
      }
    }
    for(int i = 0; i < readGroups.length; i++){
      int start = i*winStep;
      int end = start + windowLen;
      ArrayList<RefSubseq> seqs = new ArrayList<>();
      for(int j = 0; j < refAlignment.refAlns.size(); j++){
        RefSubseq refSubseq = new RefSubseq(j, refAlignment, start, end);
        if(refSubseq.start != -1){
          seqs.add(refSubseq);
        }
      }
      refSubseqs[i] = seqs;
      BufferedWriter writer = new BufferedWriter(new FileWriter(outPath + start + "_" + end + "_reads.fasta"));
      for(int j = 0; j < readGroups[i].size(); j++){
        int k = readGroups[i].getQuick(j);
        MappedRead read = mappedData.mappedReads.get(k);
        if(read.seq.length >= minLen){
          writer.write(">Read_" + k + "\n");
          writer.write(new String(read.seq) + "\n");
        }
      }
      writer.close();
    }
    logger.printf("Reads after filtering by length >= %d: %d\n", minLen, readNumAfterLenFilter);
  }

  private class UclustLauncher implements Runnable{

    private int windowID;

    UclustLauncher(int id){
      windowID = id;
    }

    @Override
    public void run(){
      try{
        int start = windowID*winStep;
        int end = start + windowLen;
        Process p = Runtime.getRuntime().exec(pathToUclust + " -threads 1 -sort length" +
            " -cluster_fast " + outPath + start + "_" + end + "_reads.fasta -id " + uclustIdentity +
            " -centroids " + outPath + "centroids_reads_" + start + "_" + end + ".fasta -uc " +
            outPath + "clusters_reads_" + start + "_" + end + ".uc");
        BufferedReader streamReader = new BufferedReader(new InputStreamReader(p.getInputStream()));
        String line;
        while((line = streamReader.readLine()) != null){
          //System.out.println(line);
        }
        streamReader = new BufferedReader(new InputStreamReader(p.getErrorStream()));
        while((line = streamReader.readLine()) != null){
          //System.out.println(line);
        }
        p.waitFor();
        if(!debug){
          File file = new File(outPath + start + "_" + end + "_reads.fasta");
          file.delete();
        }
      }catch(Exception e){
        e.printStackTrace();
        System.exit(1);
      }
    }
  }

  private void runUclust() throws InterruptedException{
    ExecutorService executor = Executors.newFixedThreadPool(threadNum);
    for(int i = 0; i < fileNum; i++) {
      UclustLauncher task = new UclustLauncher(i);
      executor.execute(task);
    }
    executor.shutdown();
    executor.awaitTermination(10, TimeUnit.HOURS);
  }

  private class VSearchLauncher implements Runnable{

    private int windowID;

    VSearchLauncher(int id){
      windowID = id;
    }

    @Override
    public void run(){
      try{
        int start = windowID*winStep;
        int end = start + windowLen;
        Process p = Runtime.getRuntime().exec(pathToUclust + " -threads 1 " +
            " -cluster_fast " + outPath + start + "_" + end + "_reads.fasta -id " + uclustIdentity +
            " -centroids " + outPath + "centroids_reads_" + start + "_" + end + ".fasta -uc " +
            outPath + "clusters_reads_" + start + "_" + end + ".uc");
        BufferedReader streamReader = new BufferedReader(new InputStreamReader(p.getInputStream()));
        String line;
        while((line = streamReader.readLine()) != null){
          //System.out.println(line);
        }
        streamReader = new BufferedReader(new InputStreamReader(p.getErrorStream()));
        while((line = streamReader.readLine()) != null){
          //System.out.println(line);
        }
        p.waitFor();
        if(!debug){
          File file = new File(outPath + start + "_" + end + "_reads.fasta");
          file.delete();
        }
      }catch(Exception e){
        e.printStackTrace();
        System.exit(1);
      }
    }
  }

  private void runVSearch() throws InterruptedException{
    ExecutorService executor = Executors.newFixedThreadPool(threadNum);
    for(int i = 0; i < fileNum; i++) {
      VSearchLauncher task = new VSearchLauncher(i);
      executor.execute(task);
    }
    executor.shutdown();
    executor.awaitTermination(10, TimeUnit.HOURS);
  }

  private void readClustersAndBuildContigs() throws IOException, InterruptedException{
    ArrayList<Cluster>[] clustersRaw = new ArrayList[fileNum];
    for(int i = 0; i < fileNum; i++){
      int start = i*winStep;
      int end = start + windowLen;
      clustersRaw[i] = new ArrayList<>();
      HashMap<String, Cluster> clusterMap = new HashMap<>();
      String line;
      BufferedReader reader = new BufferedReader(new FileReader(outPath + "clusters_reads_" +
          start + "_" + end + ".uc"));
      //TIntArrayList group = readGroups[i];
      while((line = reader.readLine()) != null){
        String[] tokens = line.split("\t");
        switch(tokens[0]){
          case "S":
            String name = tokens[8];
            int id = Integer.parseInt(tokens[1]);
            Cluster cluster = new Cluster(id, start, end, name);
            cluster.windowID = i;
            int k = Integer.parseInt(
                cluster.centroidID.substring(cluster.centroidID.indexOf('_') + 1));
            MappedRead read = mappedData.mappedReads.get(k);
            //path.reads.add(read);
            cluster.addRead(k);
            //cluster.readIDs.set(k);
            cluster.centroidSeq = new String(read.seq);
            clusterMap.put(tokens[1], cluster);
            break;
          case "H":
            Cluster cluster1 = clusterMap.get(tokens[1]);
            name = tokens[8];
            int k1 = Integer.parseInt(
                name.substring(name.indexOf('_') + 1));
            //cluster1.reads.add(mappedData.get(k));
            cluster1.addRead(k1);
            //cluster1.readIDs.set(k);
            break;
        }
      }
      reader.close();
      if(!debug){
        File file = new File(outPath + "clusters_reads_" + start + "_" + end + ".uc");
        file.delete();
        file = new File(outPath + "centroids_reads_" + start + "_" + end + ".fasta");
        file.delete();
      }
      clustersRaw[i].addAll(clusterMap.values());
    }
    ContigBuilder contigBuilder = new ContigBuilder(config, mappedData.mappedReads);
    clusters = contigBuilder.buildContigs(clustersRaw);
  }

  private void mergeContigs() throws IOException, InterruptedException{
    long time = System.currentTimeMillis();
    //graph = (Graph)graphConstructor.newInstance(config, clusters, mappedData);
    graph = new Graph(config, clusters, mappedData, windowLen, winStep, refSubseqs);
    logger.println("Graph building time: " + (System.currentTimeMillis() - time)/1000);
    time = System.currentTimeMillis();
    longestPaths = longestPathFinder.findLongestPaths(graph.vertices);
    logger.println("Longest paths time: " + (System.currentTimeMillis() - time)/1000);
    if(debug){
      debugPrinter.printLongestPaths(longestPaths);
    }
    //graph.printVertices(pathToClusters);
  }

  private class SWScoreCounter implements Runnable{

    private int vertexID;
    private ArrayList<Reference> references;
    private int[][]clusterIDToScores;


    SWScoreCounter(int[][]clusterIDToScores, ArrayList<Reference> refs, int vertexID){
      this.clusterIDToScores = clusterIDToScores;
      this.vertexID = vertexID;
      references = refs;
    }

    @Override
    public void run(){
      int refNum = references.size();
      Cluster vertex = graph.vertices[vertexID];
      byte[] seq = vertex.contig.seq;
      for(int j = 0; j < refNum; j++){
        Reference ref = references.get(j);
        int start = ref.alnToSeqStart(vertex.contig.alnStart);
        if(start == -1){
          continue;
        }
        int end = ref.alnToSeqEnd(vertex.contig.alnEnd);
        if(end == -1){
          continue;
        }
        Alignment aln = aligner.align(seq, ref.seq.substring(start, end).getBytes());
        clusterIDToScores[vertexID][j] = aln.score;
      }
    }
  }

  private class ReadScoreComputer implements Runnable{

    private int from;
    private int to;
    private int[][] readIDToScore;
    private TIntArrayList[] readIDToCluster;
    private int[][] clusterIdToScore;
    private ArrayList<Reference> refs;

    public ReadScoreComputer(ArrayList<Reference> refs, int[][] readIDToScore, TIntArrayList[] readIDToCluster, int[][] clusterIdToScore, int from, int to){
      this.refs = refs;
      this.readIDToScore = readIDToScore;
      this.readIDToCluster = readIDToCluster;
      this.clusterIdToScore = clusterIdToScore;
      this.from = from;
      this.to = to;
    }

    @Override
    public void run(){
      for(int i = from; i < to; i++){
        MappedRead read = mappedData.mappedReads.get(i);
        TIntArrayList clusterIDs = readIDToCluster[i];
        int bestClusterID = 0;
//        float maxIdentity = 0;
//        int maxReadNum = 0;
        int longestLength = 0;
        for(TIntIterator iter = clusterIDs.iterator(); iter.hasNext();){
          int clusterID = iter.next();
          Cluster cluster = graph.vertices[clusterID];
//          if(cluster.readNum > maxReadNum){
//            maxReadNum = cluster.readNum;
//            bestClusterID = clusterID;
//          }
          if(cluster.contig.seq.length > longestLength){
            longestLength = cluster.contig.seq.length;
            bestClusterID = clusterID;
          }
//          MappedRead mRead = cluster.mReads.get(i);
//          float identity = (float)mRead.aln.identity/mRead.aln.length;
//          if(identity > maxIdentity){
//            maxIdentity = identity;
//            bestClusterID = clusterID;
//          }
        }
//        if(longestLength > 0){//if(maxIdentity > identityThreshold){
//          Alignment aln = aligner.align(read.seq, graph.vertices[bestClusterID].contig.seq);
//          maxIdentity = (float)aln.identity/aln.length;
//        }
        if(longestLength > 0){//maxIdentity > identityThreshold
          readIDToScore[i] = clusterIdToScore[bestClusterID];
        }
//        else{
//          // we have to align this read;
//          for(int j = 0; j < refs.size(); j++){
//            Reference ref = refs.get(j);
//            MappedRead mRead = counter.mapReadToRegion(read.seq, ref.index, ref.contigEnds);
//            Alignment aln = aligner.align(read.seq, Arrays.copyOfRange(ref.seqB, mRead.start, mRead.end));
//            readIDToScore[i][j] = aln.score;
//          }
//        }
      }
    }

  }

  private ArrayList<Reference> sortReadsToSelectedRef(HashMap<String, int[]> selectedRef)
      throws InterruptedException{
    ArrayList<Reference> res = new ArrayList<>();
    TObjectIntHashMap<String> refNameToID = new TObjectIntHashMap<>();
    for(Map.Entry<String, int[]> entry: selectedRef.entrySet()){
      Reference seq = refAlignment.refSeqs.get(entry.getKey());
      seq.reads = new ArrayList<>();
      seq.problemReads = new ArrayList<>();
      //seq.buildIndex(config);
      res.add(seq);
      refNameToID.put(seq.name, res.size() - 1);
    }
    ExecutorService executor = Executors.newFixedThreadPool(threadNum);
    int readNum = mappedData.mappedReads.size();
    TIntArrayList[] readIDToCluster = new TIntArrayList[readNum];
    for(int i = 0; i < readIDToCluster.length; i++){
      readIDToCluster[i] = new TIntArrayList();
    }
    int refNum = selectedRef.size();
    int[][]clusterIDToScores = new int[graph.vertices.length][refNum];
    for(int i = 0; i < graph.vertices.length; i++){
      Cluster vertex = graph.vertices[i];
      for(TIntIterator iter = vertex.readIDs.iterator(); iter.hasNext();){
        readIDToCluster[iter.next()].add(i);
      }
      SWScoreCounter task = new SWScoreCounter(clusterIDToScores, res, i);
      executor.execute(task);
    }
    executor.shutdown();
    executor.awaitTermination(10, TimeUnit.HOURS);
    int[][] readIDToScore = new int[readNum][refNum];
    executor = Executors.newFixedThreadPool(threadNum);
    ArrayList<ReadScoreComputer> tasks = new ArrayList<>();
    for(int i = 0; i < readNum; i += batchSize){
      ReadScoreComputer task = new ReadScoreComputer(res, readIDToScore, readIDToCluster, clusterIDToScores, i,
          Math.min(i + batchSize, readNum));
      tasks.add(task);
      executor.execute(task);
    }
    executor.shutdown();
    executor.awaitTermination(10, TimeUnit.HOURS);
    for(int i = 0; i < mappedData.mappedReads.size(); i++){
      MappedRead mappedRead = mappedData.mappedReads.get(i);
      int pairID = mappedData.pairs[i];
      PairedRead pairedRead;
      TIntArrayList refIDs = new TIntArrayList();
      if(pairID == -1){
        byte[] s;
        if(mappedRead.reverse == 1){
          s = getComplement(mappedRead.seq);
        }else{
          s = mappedRead.seq;
        }
        if(mappedRead.n == 1){
          pairedRead = new PairedRead(mappedRead.name, s, NULL_SEQ, mappedRead.q, NULL_SEQ);
        }else{
          pairedRead = new PairedRead(mappedRead.name, NULL_SEQ, s, NULL_SEQ, mappedRead.q);
        }
        for(Reference ref: res){
          ref.problemReads.add(new PairedRead(pairedRead));
        }
////        if(readIDToCluster[i].size() == 0){//mappedRead.seq.length < minLen ||
////          continue;
////        }
//        int maxScore = 0;
//        for(int j = 0; j < refNum; j++){
//          int sumScore = readIDToScore[i][j];
//          if(sumScore > maxScore){
//            maxScore = sumScore;
//            refIDs.clear();
//            refIDs.add(j);
//          }else if(sumScore == maxScore){
//            refIDs.add(j);
//          }
//        }
      }else{
        if(pairID < i){
          continue;
        }
        MappedRead readPair = mappedData.mappedReads.get(pairID);
        byte[] s1;
        if(mappedRead.reverse == 1){
          s1 = getComplement(mappedRead.seq);
        }else{
          s1 = mappedRead.seq;
        }
        byte[] s2;
        if(readPair.reverse == 1){
          s2 = getComplement(readPair.seq);
        }else{
          s2 = readPair.seq;
        }
        if(mappedRead.n == 1){
          pairedRead = new PairedRead(mappedRead.name, s1, s2, mappedRead.q, readPair.q);
        }else{
          pairedRead = new PairedRead(mappedRead.name, s2, s1, readPair.q, mappedRead.q);
        }
        if(mappedData.isConcordant[i] == 0 || readIDToCluster[i].size() == 0 || readIDToCluster[pairID].size() == 0){//mappedRead.seq.length < minLen || readPair.seq.length < minLen ||
          for(Reference ref: res){
            ref.problemReads.add(new PairedRead(pairedRead));
          }
          continue;
        }
        for(int j = 0, maxScore = 0; j < refNum; j++){
          int sumScore = readIDToScore[i][j] + readIDToScore[pairID][j];
          if(sumScore > maxScore){
            maxScore = sumScore;
            refIDs.clear();
            refIDs.add(j);
          }else if(sumScore == maxScore){
            refIDs.add(j);
          }
        }
        if(res.size() > 0){
          res.get(refIDs.getQuick(0)).reads.add(pairedRead);
          if(refIDs.size() > 1){
            for(int j = 1; j < refIDs.size(); j++){
              res.get(refIDs.getQuick(j)).reads.add(new PairedRead(pairedRead));
            }
          }
        }
      }

    }

    return res;
  }


  public ArrayList<Reference> selectReferences(ArrayList<PairedRead> reads) throws InterruptedException, IOException{
    long time = System.currentTimeMillis();
    mappedData = mapper.mapAllReads(reads, refAlignment);
    logger.println("Time: " + (System.currentTimeMillis() - time)/1000);
    time = System.currentTimeMillis();
    int maxReadLen = 0;
    for(MappedRead read: mappedData.mappedReads){
      int len = read.end - read.start;
      if(len > maxReadLen){
        maxReadLen = len;
      }
    }
    windowLen = 3*maxReadLen/2 + 1;
    winStep = maxReadLen/2 - 1;
    fileNum = (refAlignment.length - windowLen)/winStep + 1;
    groupReads();
    if(debug){
      logger.printTotalReadNumInGroups(readGroups);
    }
    logger.println("Preprocess time: " + (System.currentTimeMillis() - time)/1000);
    time = System.currentTimeMillis();
//    runUclust();
    runVSearch();
    logger.println("UClust time: " + (System.currentTimeMillis() - time)/1000);
    time = System.currentTimeMillis();
    readClustersAndBuildContigs();
    logger.println("Contig building time: " + (System.currentTimeMillis() - time)/1000);
    //logger.printClustersToFiles(clusters, mappedData.mappedReads);
    if(debug){
      logger.printTotalReadNumInClusters(clusters);
    }
    time = System.currentTimeMillis();
    mergeContigs();
    logger.println("Contig merging time: " + (System.currentTimeMillis() - time)/1000);
    if(debug){
      logger.printTotalReadNumInGraph(graph.vertices);
      logger.printTotalReadNumInPaths(longestPaths);
      logger.printTotalReadNumInLongPaths(longestPaths);
    }
    time = System.currentTimeMillis();
    HashMap<String, int[]> selectedRef = minRefSetSelector.chooseReferences(longestPaths);
    logger.println("Reference selection time: " + (System.currentTimeMillis() - time)/1000);
    time = System.currentTimeMillis();
    //ArrayList<Reference> res = sortReadsToSelectedRef(selectedRef, reads);
    ArrayList<Reference> res = sortReadsToSelectedRef(selectedRef);
    logger.println("Read assigning time: " + (System.currentTimeMillis() - time)/1000);
    return res;
  }

  static void addParameters(ArgumentParser parser){
    parser.description("This tool will " +
        "choose one or more references from given paired reads using given set of aligned reference " +
        "sequences. All the parameters should be provided in the config file.");
    parser.addArgument("-c").dest("config_file").help("Path to config file").required(true);
  }

  static void run(Namespace parsedArgs){
    try{
      SAXBuilder jdomBuilder = new SAXBuilder();
      Document jdomDocument = jdomBuilder.build(parsedArgs.getString("config_file"));
      ReferenceFinder refFinder = new ReferenceFinder(jdomDocument);
      DataReader dataReader = DataReader.getInstance(jdomDocument);
      refFinder.selectReferences(dataReader.pairedReads);
    }catch(Exception e){
      e.printStackTrace();
    }
  }

  public static void main(String[] args){
    ArgumentParser parser = ArgumentParsers.newFor("ReferenceFinder").build();
    addParameters(parser);
    if(args.length == 0){
      parser.printUsage();
      return;
    }
    try{
      Namespace parsedArgs = parser.parseArgs(args);
      run(parsedArgs);
    }catch(ArgumentParserException e){
      parser.handleError(e);
    }

  }

}
