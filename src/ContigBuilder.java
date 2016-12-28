import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TIntObjectHashMap;
import org.jdom2.Document;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Created with IntelliJ IDEA.
 * Date: 08.08.15
 */
class ContigBuilder extends Constants{

  KMerCounter counter;
  Aligner aligner;
  ArrayList<MappedRead> mappedReads;
  private ReferenceAlignment refAlignment;
  private int K;
  private int threadNum;

  ContigBuilder(Document document, ArrayList<MappedRead> mappedReads){
    this.mappedReads = mappedReads;
    counter = KMerCounter.getInstance(document);
    aligner = new SmithWatermanGotoh(document);
    K = counter.K;
    refAlignment = ReferenceAlignment.getInstance(document);
    threadNum = Integer.parseInt(document.getRootElement().getChildText("ThreadNumber"));
  }

  private void buildContigAln(Cluster cluster){
    int[] seqToAln = refAlignment.refSeqs.get(cluster.centroidID).seqToAln;
    TByteArrayList list = new TByteArrayList();
    TIntArrayList cov = new TIntArrayList();
    int prev = seqToAln[cluster.contig.centroidStart + cluster.centroidStart];
    list.add(cluster.contig.seq[0]);
    cov.add(cluster.contig.coverage[0]);
    for(int i = cluster.contig.centroidStart + cluster.centroidStart + 1, j = 1;
        i < cluster.contig.centroidEnd + cluster.centroidStart; i++, j++){
      int k = seqToAln[i];
      for(int n = 0; n < k - prev - 1; n++){
        list.add(GAP);
        cov.add(0);
      }
      list.add(cluster.contig.seq[j]);
      cov.add(cluster.contig.coverage[j]);
      prev = k;
    }
    cluster.contig.seqAln = list.toArray();
    cluster.contig.coverageAln = cov.toArray();
//    if(cluster.contig.seqAln.length != cluster.contig.alnEnd - cluster.contig.alnStart){
//      int a = 0;
//    }
  }

  private class ContigBuilding implements Runnable{

    private Cluster cluster;
    public LinkedList<Cluster> res;

    public ContigBuilding(Cluster cluster){
      this.cluster = cluster;
    }

    private int[][] countReads(){
      int centroidLen = cluster.centroidSeq.length();
      int[][] counts = new int[centroidLen][5];
      TIntIterator iter = cluster.readIDs.iterator();
      while(iter.hasNext()){
        int id = iter.next();
        //for(int id = readIDs.nextSetBit(0); id >= 0; id = readIDs.nextSetBit(id+1)){
        MappedRead read = mappedReads.get(id);
        MappedRead mRead = counter.mapReadToRegion(read.seq, cluster.centroidIndex, 0, centroidLen);
        if(mRead == null){
          continue;
        }
        byte[] s2 = new byte[mRead.end - mRead.start];
        for(int i = 0; i < s2.length; i++){
          s2[i] = (byte)cluster.centroidSeq.charAt(i + mRead.start);
        }
        Alignment alignment = aligner.align(read.seq, s2);
        mRead.aln = alignment;
        //mRead.name = read.name;
        byte[] seq1 = alignment.sequence1;
        int k = mRead.start + alignment.start2;
        for(int j = 0; j < alignment.length; j++){
          if(seq1[j] < (byte)'a'){
            counts[k][nToI[seq1[j]]] ++;
            k ++;
          }
        }
        cluster.mReads.put(id, mRead);
      }
      return counts;
    }

    @Override
    public void run(){

      boolean centroidIsRef = false;
      int[] seqToAln = null;
      if(cluster.refIDs.size() > 0){
        // centroid is reference -> compute contig start and end
        centroidIsRef = true;
        Reference seq = refAlignment.refSeqs.get(cluster.centroidID);
        if(seq == null){
          RefSubseq refSubseq = cluster.refSubseqs.get(0);
          seq = refAlignment.refAlns.get(refSubseq.refID);
          cluster.centroidID = seq.name;
          cluster.centroidStart = refSubseq.start;
          cluster.centroidEnd = refSubseq.end;
          cluster.centroidSeq = seq.seq.substring(cluster.centroidStart, cluster.centroidEnd);
        }
        seqToAln = seq.seqToAln;
      }
      int centroidLen = cluster.centroidSeq.length();
      Contig contig = new Contig();
      cluster.contig = contig;
      cluster.mReads = new TIntObjectHashMap<>();
      cluster.buildCentroidIndex(K);
      int[][] counts = countReads();
      res = new LinkedList<>();
      TByteArrayList builder = new TByteArrayList();
      TIntArrayList cover = new TIntArrayList();
      boolean started = false;
      int contigNum = 0;
      for(int j = 0; j < centroidLen; j++){
        int maxCount = counts[j][4];
        int maxCountIndex = 4;
        for(int k = 0; k < 4; k++){
          if(counts[j][k] > maxCount){
            maxCount = counts[j][k];
            maxCountIndex = k;
          }
        }
        if(started){
          if(maxCount == 0){
            if(centroidIsRef){
              contig.alnEnd = seqToAln[cluster.centroidStart - 1 + j] + 1;
            }
            started = false;
            contig.centroidEnd = j;
            Cluster smallCluster = new Cluster(cluster, builder.toArray(), cover.toArray());
            smallCluster.name += "_" + Integer.toString(contigNum);
            contigNum ++;
            if(smallCluster.refIDs.size() > 0){
              buildContigAln(smallCluster);
            }
            res.add(smallCluster);
            builder = new TByteArrayList();
            cover = new TIntArrayList();
          }else{
            builder.add(iToN[maxCountIndex]);
            cover.add(maxCount);
          }
        }else{
          if(maxCount > 0){
            started = true;
            if(centroidIsRef){
              contig.alnStart = seqToAln[cluster.centroidStart + j];
            }
            contig.centroidStart = j;
            builder.add(iToN[maxCountIndex]);
            cover.add(maxCount);
          }
        }
      }
      if(started){
        if(centroidIsRef){
          contig.alnEnd = seqToAln[cluster.centroidEnd - 1] + 1;
        }
        contig.centroidEnd = centroidLen;
        Cluster smallCluster = new Cluster(cluster, builder.toArray(), cover.toArray());
        smallCluster.name += "_" + Integer.toString(contigNum);
        contigNum ++;
        if(smallCluster.refIDs.size() > 0){
          buildContigAln(smallCluster);
        }
        res.add(smallCluster);
      }
      if(res.size() == 0){
//        contig.centroidStart = 0;
//        contig.centroidEnd = cluster.centroidSeq.length();
//        contig.coverage = new int[cluster.centroidSeq.length()];
//        contig.seqB = cluster.centroidSeq.getBytes();
        cluster.contig = null;
        res.add(cluster);
      }
    }

  }

  private class ContigReBuilding implements Runnable{

    private Cluster cluster;

    public ContigReBuilding(Cluster cluster){
      this.cluster = cluster;
    }

    private int[][] countReads(){
      int centroidLen = cluster.centroidSeq.length();
      int[][] counts = new int[centroidLen][5];
      TIntIterator iter = cluster.readIDs.iterator();
      while(iter.hasNext()){
        int id = iter.next();
        //for(int id = readIDs.nextSetBit(0); id >= 0; id = readIDs.nextSetBit(id+1)){
        MappedRead mRead;
        if(cluster.readsToRemap.contains(id)){
          MappedRead read = mappedReads.get(id);
          mRead = counter.mapReadToRegion(read.seq, cluster.centroidIndex, 0, centroidLen);
          if(mRead == null){
            continue;
          }
          byte[] s2 = new byte[mRead.end - mRead.start];
          for(int i = 0; i < s2.length; i++){
            s2[i] = (byte) cluster.centroidSeq.charAt(i + mRead.start);
          }
          Alignment alignment = aligner.align(read.seq, s2);
          mRead.aln = alignment;
          cluster.mReads.put(id, mRead);
        }else{
          mRead = cluster.mReads.get(id);
        }
        //mRead.name = read.name;
        byte[] seq1 = mRead.aln.sequence1;
        int k = mRead.start + mRead.aln.start2;
        for(int j = 0; j < mRead.aln.length; j++){
          if(seq1[j] < (byte)'a'){
            counts[k][nToI[seq1[j]]] ++;
            k ++;
          }
        }
      }
      return counts;
    }

    @Override
    public void run(){

      Reference seq = refAlignment.refSeqs.get(cluster.centroidID);
      if(seq == null){
        RefSubseq refSubseq = cluster.refSubseqs.get(0);
        seq = refAlignment.refAlns.get(refSubseq.refID);
        cluster.centroidID = seq.name;
        cluster.centroidStart = refSubseq.start;
        cluster.centroidEnd = refSubseq.end;
        cluster.centroidSeq = seq.seq.substring(cluster.centroidStart, cluster.centroidEnd);
      }
      int[] seqToAln = seq.seqToAln;
      int centroidLen = cluster.centroidSeq.length();
      Contig contig = new Contig();
      cluster.contig = contig;
      cluster.buildCentroidIndex(K);
      int[][] counts = countReads();
      TByteArrayList builder = new TByteArrayList();
      TIntArrayList cover = new TIntArrayList();
      boolean started = false;
      for(int j = 0; j < centroidLen; j++){
        int maxCount = counts[j][4];
        int maxCountIndex = 4;
        for(int k = 0; k < 4; k++){
          if(counts[j][k] > maxCount){
            maxCount = counts[j][k];
            maxCountIndex = k;
          }
        }
        if(started){
          if(maxCount == 0){
            contig.alnEnd = seqToAln[cluster.centroidStart - 1 + j] + 1;
            started = false;
            contig.centroidEnd = j;
            break;
          }else{
            builder.add(iToN[maxCountIndex]);
            cover.add(maxCount);
          }
        }else{
          if(maxCount > 0){
            started = true;
            contig.alnStart = seqToAln[cluster.centroidStart + j];
            contig.centroidStart = j;
            builder.add(iToN[maxCountIndex]);
            cover.add(maxCount);
          }
        }
      }
      if(started){
        contig.alnEnd = seqToAln[cluster.centroidEnd - 1] + 1;
        contig.centroidEnd = centroidLen;
      }
      contig.seq = builder.toArray();
      contig.coverage = cover.toArray();
      buildContigAln(cluster);
    }

  }

  LinkedList<Cluster>[] buildContigs(ArrayList<Cluster>[] clusters)
      throws InterruptedException {
    ArrayList<ContigBuilding>[] tasks = new ArrayList[clusters.length];
    ExecutorService executor = Executors.newFixedThreadPool(threadNum);
    for(int i = 0; i < clusters.length; i++) {
      tasks[i] = new ArrayList<>();
      for (Cluster cluster : clusters[i]) {
        ContigBuilding task = new ContigBuilding(cluster);
        executor.execute(task);
        tasks[i].add(task);
      }
    }
    executor.shutdown();
    executor.awaitTermination(10, TimeUnit.HOURS);
    LinkedList<Cluster>[] res = new LinkedList[clusters.length];
    for(int i = 0; i < clusters.length; i++) {
      res[i] = new LinkedList<>();
      for (ContigBuilding task : tasks[i]) {
        res[i].addAll(task.res);
      }
    }
    return res;
  }

  LinkedList<LinkedList<Cluster>> buildContigs(Collection<Cluster> clusters)
      throws InterruptedException {
    ArrayList<ContigBuilding> tasks = new ArrayList<>();
    ExecutorService executor = Executors.newFixedThreadPool(threadNum);
    for(Cluster cluster: clusters) {
      ContigBuilding task = new ContigBuilding(cluster);
      executor.execute(task);
      tasks.add(task);
    }
    executor.shutdown();
    executor.awaitTermination(10, TimeUnit.HOURS);
    LinkedList<LinkedList<Cluster>> res = new LinkedList<>();
    for(ContigBuilding task : tasks) {
      res.add(task.res);
    }
    return res;
  }

  void rebuildContigs(Collection<Cluster> clusters)
      throws InterruptedException {
    ExecutorService executor = Executors.newFixedThreadPool(threadNum);
    for(Cluster cluster: clusters) {
      executor.execute(new ContigReBuilding(cluster));
    }
    executor.shutdown();
    executor.awaitTermination(10, TimeUnit.HOURS);
    for(Cluster cluster: clusters){
      cluster.readsToRemap.clear();
    }
  }


}
