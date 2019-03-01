import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.iterator.TObjectIntIterator;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TIntFloatHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import org.jdom2.Document;
import org.jdom2.Element;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Created by Gennady on 08.09.2015.
 */
class Graph extends Constants{

  private int minReadNum;
  private float uclustIdentity;

  Cluster[] vertices;
  private KMerCounter counter;
  private Aligner aligner;
  private ReferenceAlignment refAlignment;
  private MappedData mappedData;
  private Document config;
  private Logger logger;
  private EdgeCreator edgeCreator;
  private float vertexWeight;
  private int[][] readIDToRefIDs;
  private String outPath;

  private ArrayList<RefSubseq>[] refSubseqs;
  private float delta;
  private float similarityThreshold;
  private int threadNum;
  private boolean debug;

  private float computeIdentity(Cluster cluster_i, Cluster cluster_j){
    int identity = 0;
    int len = 0;
    for(int i = cluster_j.contig.alnStart - cluster_i.contig.alnStart, j = 0;
        i < cluster_i.contig.seqAln.length && j < cluster_j.contig.seqAln.length; i++, j++){
      if(cluster_i.contig.seqAln[i] == GAP){
        if(cluster_j.contig.seqAln[j] != GAP){
          len++;
        }
      }else{
        if(cluster_j.contig.seqAln[j] == GAP){
          len++;
        }else{
          if(cluster_i.contig.seqAln[i] == cluster_j.contig.seqAln[j]){
            identity++;
          }
          len++;
        }
      }
    }
    return (float) identity/len;
  }

  private void trimUncoveredEnds(){
    for(Cluster cluster : vertices){
      int startNew = 0;
      int endNew = cluster.contig.seq.length;
      if(cluster.neighborsBackward.size() == 0){
        // starts path
        if(cluster.neighborsForward.size() == 0){
          cluster.markedForDeletion = true;
          continue;
        }
        for(int i = 0; i < cluster.contig.coverage.length; i++){
          if(cluster.contig.coverage[i] < minReadNum){
            startNew++;
          }else{
            break;
          }
        }
        if(startNew == endNew){
          cluster.markedForDeletion = true;
          continue;
        }
        if(startNew > 0){
          cluster.contig.coverage =
              Arrays.copyOfRange(cluster.contig.coverage, startNew, endNew);
          cluster.contig.seq =
              Arrays.copyOfRange(cluster.contig.seq, startNew, endNew);
          int[] seqToAln = refAlignment.refSeqs.get(cluster.centroidID).seqToAln;
          cluster.contig.centroidStart += startNew;
          cluster.contig.centroidEnd = cluster.contig.centroidStart + endNew - startNew;
          int newAlnStart = seqToAln[cluster.centroidStart + cluster.contig.centroidStart];
          int newAlnEnd = seqToAln[cluster.centroidStart + cluster.contig.centroidEnd - 1] + 1;
          cluster.contig.seqAln = Arrays.copyOfRange(cluster.contig.seqAln,
              newAlnStart - cluster.contig.alnStart, newAlnEnd - cluster.contig.alnStart);
          cluster.contig.coverageAln = Arrays.copyOfRange(cluster.contig.coverageAln,
              newAlnStart - cluster.contig.alnStart, newAlnEnd - cluster.contig.alnStart);
          cluster.contig.alnStart = newAlnStart;
          cluster.contig.alnEnd = newAlnEnd;
          ArrayList<Cluster> neighborsForward = new ArrayList<>();
          TFloatArrayList weightsForward = new TFloatArrayList();
          for(int i = 0; i < cluster.neighborsForward.size(); i++){
            Cluster neighbor = cluster.neighborsForward.get(i);
            float neighborWeight = cluster.weightsForward.getQuick(i);
            if(neighbor.contig.alnStart < cluster.contig.alnStart){
              int j = neighbor.neighborsBackward.indexOf(cluster);
              neighbor.neighborsBackward.remove(j);
              neighbor.weightsBackward.removeAt(j);
              if(cluster.contig.alnStart < neighbor.contig.alnEnd){
                neighbor.neighborsForward.add(cluster);
                neighbor.weightsForward.add(neighborWeight);
                cluster.neighborsBackward.add(neighbor);
                cluster.weightsBackward.add(neighborWeight);
              }
            }else{
              neighborsForward.add(neighbor);
              weightsForward.add(cluster.weightsForward.getQuick(i));
            }
          }
          cluster.neighborsForward = neighborsForward;
          cluster.weightsForward = weightsForward;
        }
      }else if(cluster.neighborsForward.size() == 0){
        // ends path
        for(int i = cluster.contig.coverage.length - 1; i >= 0; i--){
          if(cluster.contig.coverage[i] < minReadNum){
            endNew--;
          }else{
            break;
          }
        }
        if(endNew == 0){
          cluster.markedForDeletion = true;
          continue;
        }
        if(endNew < cluster.contig.seq.length){
          cluster.contig.coverage =
              Arrays.copyOfRange(cluster.contig.coverage, startNew, endNew);
          cluster.contig.seq =
              Arrays.copyOfRange(cluster.contig.seq, startNew, endNew);
          int[] seqToAln = refAlignment.refSeqs.get(cluster.centroidID).seqToAln;
          cluster.contig.centroidStart += startNew;
          cluster.contig.centroidEnd = cluster.contig.centroidStart + endNew - startNew;
          int newAlnStart = seqToAln[cluster.centroidStart + cluster.contig.centroidStart];
          int newAlnEnd = seqToAln[cluster.centroidStart + cluster.contig.centroidEnd - 1] + 1;
          cluster.contig.seqAln = Arrays.copyOfRange(cluster.contig.seqAln,
              newAlnStart - cluster.contig.alnStart, newAlnEnd - cluster.contig.alnStart);
          cluster.contig.coverageAln = Arrays.copyOfRange(cluster.contig.coverageAln,
              newAlnStart - cluster.contig.alnStart, newAlnEnd - cluster.contig.alnStart);
          cluster.contig.alnStart = newAlnStart;
          cluster.contig.alnEnd = newAlnEnd;
        }
      }
    }
    Arrays.sort(vertices);
  }

  private ArrayList<Cluster> mergeIncludedClusters(ArrayList<Cluster> clustersWithReads){
    for(int i = 0; i < clustersWithReads.size(); i++){
      Cluster cluster_i = clustersWithReads.get(i);
      if(cluster_i.markedForDeletion){
        continue;
      }
      for(int j = i + 1; j < clustersWithReads.size(); j++){
        Cluster cluster_j = clustersWithReads.get(j);
        if(cluster_j.contig.alnStart >= cluster_i.contig.alnEnd){
          break;
        }
        if(cluster_j.markedForDeletion){
          continue;
        }
        if(cluster_i.contig.alnEnd >= cluster_j.contig.alnEnd){
          float identity = computeIdentity(cluster_i, cluster_j);
          if(identity > uclustIdentity){
            for(TIntIterator iter = cluster_j.readIDs.iterator(); iter.hasNext(); ){
              int id = iter.next();
              if(cluster_i.readIDs.add(id)){
                cluster_i.readsToRemap.add(id);
              }
            }
            cluster_i.readNum = cluster_i.readIDs.size();
            cluster_j.markedForDeletion = true;
          }
        }
      }
    }

    ArrayList<Cluster> clustersToRebuild = new ArrayList<>();
    ArrayList<Cluster> filteredVertices = new ArrayList<>();
    for(Cluster cluster : clustersWithReads){
      if(!cluster.markedForDeletion){
        filteredVertices.add(cluster);
        if(cluster.readsToRemap.size() > 0){
          clustersToRebuild.add(cluster);
        }
      }
    }
    ContigBuilder contigBuilder = new ContigBuilder(config, mappedData.mappedReads);
    try{
      contigBuilder.rebuildContigs(clustersToRebuild);
    }catch(InterruptedException e){
      e.printStackTrace();
      System.exit(1);
    }
    Collections.sort(filteredVertices);
    return filteredVertices;
  }

  private void mergeIncludedClusters(){
    for(int i = 0; i < vertices.length; i++){
      Cluster cluster_i = vertices[i];
      if(cluster_i.markedForDeletion){
        continue;
      }
      for(int j = i + 1; j < vertices.length; j++){
        Cluster cluster_j = vertices[j];
        if(cluster_j.contig.alnStart >= cluster_i.contig.alnEnd){
          break;
        }
        if(cluster_j.markedForDeletion){
          continue;
        }
        if(cluster_i.contig.alnEnd >= cluster_j.contig.alnEnd){
          float identity = computeIdentity(cluster_i, cluster_j);
          if(identity > uclustIdentity){
            for(TIntIterator iter = cluster_j.readIDs.iterator(); iter.hasNext(); ){
              int id = iter.next();
              if(cluster_i.readIDs.add(id)){
                cluster_i.readsToRemap.add(id);
              }
            }
            cluster_i.readNum = cluster_i.readIDs.size();
            cluster_j.markedForDeletion = true;
          }
        }
      }
    }

    ArrayList<Cluster> clustersToRebuild = new ArrayList<>();
    ArrayList<Cluster> filteredVertices = new ArrayList<>();
    for(Cluster cluster : vertices){
      if(!cluster.markedForDeletion){
        filteredVertices.add(cluster);
        if(cluster.readsToRemap.size() > 0){
          cluster.edgesNeedToBeRepaired = true;
          clustersToRebuild.add(cluster);
        }
      }else{
        for(Cluster backNeighbor: cluster.neighborsBackward){
          int i = backNeighbor.neighborsForward.indexOf(cluster);
          backNeighbor.neighborsForward.remove(i);
          backNeighbor.weightsForward.removeAt(i);
        }
        for(Cluster forwardNeighbor: cluster.neighborsForward){
          int i = forwardNeighbor.neighborsBackward.indexOf(cluster);
          forwardNeighbor.neighborsBackward.remove(i);
          forwardNeighbor.weightsBackward.removeAt(i);
        }
      }
    }
    ContigBuilder contigBuilder = new ContigBuilder(config, mappedData.mappedReads);
    try{
      contigBuilder.rebuildContigs(clustersToRebuild);
    }catch(InterruptedException e){
      e.printStackTrace();
      System.exit(1);
    }
    vertices = filteredVertices.toArray(new Cluster[filteredVertices.size()]);
    Arrays.sort(vertices);
    edgeCreator.repairEdges(vertices);
  }

  Graph(Document document, LinkedList<Cluster>[] clusters, MappedData mappedData,
        int windowLen, int winStep, ArrayList<RefSubseq>[] refSubseqs){
    config = document;
    logger = Logger.getInstance(document);
    this.refSubseqs = refSubseqs;
    Element graphElement = document.getRootElement().getChild("ReferenceSelector").getChild("Graph");
    debug = Boolean.parseBoolean(graphElement.getChildText("Debug"));
    delta = Float.parseFloat(document.getRootElement().getChild("ReferenceSelector").getChildText("Delta"));
    minReadNum = Integer.parseInt(graphElement.getChildText("MinReadNumber"));
    vertexWeight = Float.parseFloat(graphElement.getChildText("VertexWeight"));
    similarityThreshold = Float.parseFloat(graphElement.getChildText("SimilarityThreshold"));
    uclustIdentity = Float.parseFloat(
        document.getRootElement().getChild("ReferenceSelector").getChildText("UclustIdentity"));
    refAlignment = ReferenceAlignment.getInstance(document);
    this.mappedData = mappedData;
    outPath = document.getRootElement().getChildText("OutPath");
    edgeCreator = new CommonReferenceAndReadEdgeCreator();
    aligner = new SmithWatermanGotoh(document);
    counter = KMerCounter.getInstance(document);
    threadNum = Integer.parseInt(document.getRootElement().getChildText("ThreadNumber"));
    if(threadNum == -1){
      threadNum = Runtime.getRuntime().availableProcessors();
    }
    // merge vertices with reads as centroids to some vertices with reference centroids
    assignReferencesToClusters(clusters);

    if(debug){
      logger.println("Read only cluster merging done!");
      logger.printClusterStats(clusters, windowLen, winStep);
    }

    ArrayList<Cluster> clustersWithReadsAndReferences = new ArrayList<>();
    for(LinkedList<Cluster> list: clusters){
      for(Cluster cluster: list){
        if(cluster.refNum > 0){
          clustersWithReadsAndReferences.add(cluster);
        }
      }
    }
    Collections.sort(clustersWithReadsAndReferences);
    ArrayList<Cluster> filteredClusters = mergeIncludedClusters(clustersWithReadsAndReferences);

    logger.println("Included cluster treatment done!");
    if(debug){
      logger.printTotalReadNumInClusters(filteredClusters);
      logger.printClusterStats(filteredClusters, windowLen, winStep, clusters.length);
    }

    vertices = filteredClusters.toArray(new Cluster[filteredClusters.size()]);
    edgeCreator.createEdges(this.vertices);

    logger.println("Edges created!");
    trimUncoveredEnds();
    mergeIncludedClusters();
    trimUncoveredEnds();
    logger.println("Low covered ends trimming done!");
    if(debug){
      logger.printTotalReadNumInGraph(vertices);
      logger.printVertexStats(vertices, windowLen, winStep, clusters.length);
    }

    assignClusterWeights();
  }

  private void assignClusterWeights(){
    ArrayList<Cluster> largeClusters = new ArrayList<>();
    for(Cluster cluster : vertices){
      if(!cluster.markedForDeletion){
        largeClusters.add(cluster);
        if(cluster.readNum >= minReadNum){
          cluster.weight = vertexWeight;
        }
      }else{
        for(Cluster frontNeighbor : cluster.neighborsForward){
          if(!frontNeighbor.markedForDeletion){
            int i = frontNeighbor.neighborsBackward.indexOf(cluster);
            frontNeighbor.neighborsBackward.remove(i);
            frontNeighbor.weightsBackward.removeAt(i);
          }
        }
        for(Cluster backNeighbor : cluster.neighborsBackward){
          if(!backNeighbor.markedForDeletion){
            int i = backNeighbor.neighborsForward.indexOf(cluster);
            backNeighbor.neighborsForward.remove(i);
            backNeighbor.weightsForward.removeAt(i);
          }
        }
      }
    }
    vertices = largeClusters.toArray(new Cluster[largeClusters.size()]);
  }

  private class PairScoring implements Runnable, Comparable{

    public Cluster cluster;
    public int[] refCounts;
    public float[] countsNorm;

    public PairScoring(Cluster cluster){
      this.cluster = cluster;
    }

    @Override
    public void run(){
      refCounts = new int[refAlignment.refNum];
      countsNorm = new float[refAlignment.refNum];
      for(TIntIterator iter = cluster.readIDs.iterator(); iter.hasNext();){
        int pairID = mappedData.pairs[iter.next()];
        if(pairID == -1){
          continue;
        }
        int[] refIDs = readIDToRefIDs[pairID];
        if(refIDs != null){
          for(int refID : refIDs){
            refCounts[refID]++;
          }
        }
      }
      float sumSqr = 0;
      for(int refCount: refCounts){
        sumSqr += refCount*refCount;
      }
      sumSqr = (float)Math.sqrt(sumSqr);
      for(int i = 0; i < refCounts.length; i++){
        countsNorm[i] = (float)refCounts[i]/sumSqr;
      }
    }

    @Override
    public int compareTo(Object o){
      PairScoring scoring = (PairScoring) o;
      return scoring.cluster.readNum - cluster.readNum;
    }
  }

  private PairScoring[][] computePairCounts(LinkedList<Cluster>[] clusters){
    readIDToRefIDs = new int[mappedData.mappedReads.size()][];
    for(LinkedList<Cluster> clusterList: clusters){
      for(Cluster cluster: clusterList){
        int[] refIDs = cluster.refIDs.toArray();
        for(TIntIterator iter = cluster.readIDs.iterator(); iter.hasNext(); ){
          readIDToRefIDs[iter.next()] = refIDs;
        }
      }
    }
    PairScoring[][] tasks = new PairScoring[clusters.length][];
    try{
      ExecutorService executor = Executors.newFixedThreadPool(threadNum);
      for(int i = 0; i < clusters.length; i++){
        LinkedList<Cluster> clusterList = clusters[i];
        tasks[i] = new PairScoring[clusterList.size()];
        int j = 0;
        for(Cluster cluster: clusterList){
          PairScoring task = new PairScoring(cluster);
          tasks[i][j++] = task;
          if(!cluster.markedForDeletion){
            executor.execute(task);
          }
        }
      }
      executor.shutdown();
      executor.awaitTermination(10, TimeUnit.HOURS);

      if(debug){
        for(PairScoring[] tasks_i : tasks){
          for(PairScoring task : tasks_i){
            if(task.cluster.markedForDeletion){
              continue;
            }
            BufferedWriter writer = new BufferedWriter(new FileWriter(outPath +
                task.cluster.name + ".cluster"));
            StringBuilder builder = new StringBuilder();
            builder.append(task.cluster.name).append('\n');
            builder.append("Total reads:").append(task.cluster.readNum).append('\n');

            for(TIntIterator iter = task.cluster.readIDs.iterator(); iter.hasNext(); ){
              int id = iter.next();
              builder.append(id).append('\n');
            }
            builder.append("Total references:").append(task.cluster.refNum).append('\n');
            for(TIntIterator iter = task.cluster.refIDs.iterator(); iter.hasNext(); ){
              Reference s = refAlignment.refAlns.get(iter.next());
              builder.append(s.name).append('\n');
            }
            builder.append("Ref counts:\n");
            for(int i = 0; i < task.refCounts.length; i++){
              int refCount = task.refCounts[i];
              String name = refAlignment.refAlns.get(i).name;
              builder.append(name).append('\t').append(refCount).append('\n');
            }
            builder.append("Contigs stats: \n");
            builder.append(task.cluster.name).append('\n');
            builder.append("Read num: ").append(task.cluster.readNum).append('\n');
            builder.append("Contig len = ").append(task.cluster.contig.seq.length).append('\n');
            float avIdentity = 0;
            for(TIntObjectIterator<MappedRead> iter = task.cluster.mReads.iterator(); iter.hasNext(); ){
              iter.advance();
              MappedRead mRead = iter.value();
              if(mRead.aln != null){
                avIdentity += (float) mRead.aln.identity/mRead.aln.length;
              }
            }
            avIdentity /= task.cluster.readNum;
            builder.append(" Av identity = ").append(avIdentity);
            writer.write(builder.toString());
            writer.close();
          }
        }
      }
    }catch(Exception e){
      e.printStackTrace();
      System.exit(1);
    }
    return tasks;
  }

  private class SWScore implements Runnable{

    Cluster cluster;
    float score;
    float identity;
    float coverage;
    private boolean changeCentroid;

    SWScore(Cluster cluster, boolean changeCentroid){
      this.cluster = cluster;
      this.changeCentroid = changeCentroid;
    }

    private class RefScore implements Comparable{
      public RefSubseq ref;
      public float score;
      public float identity;

      public RefScore(RefSubseq refSubseq, float s, float id, float cov){
        ref = refSubseq;
        score = s;
        identity = id;
        coverage = cov;
      }

      @Override
      public int compareTo(Object o){
        RefScore s = (RefScore)o;
        if(score == s.score){
          return 0;
        }
        return (score > s.score)?-1:1;
      }
    }

    @Override
    public void run(){
      byte[] seq = new byte[cluster.contig.seq.length];
      int len = 0;
      for(int i = 0; i < cluster.contig.seq.length; i++){
        byte ch = cluster.contig.seq[i];
        if(ch != GAP){
          seq[len] = ch;
          len++;
        }
      }
      seq = Arrays.copyOfRange(seq, 0, len);
      int fullScore = len*aligner.match;
      ArrayList<RefScore> scores = new ArrayList<>();
      for(RefSubseq ref: refSubseqs[cluster.windowID]){
        Reference reference = refAlignment.refAlns.get(ref.refID);
        String refSeq = reference.seq.substring(ref.start, ref.end);
        try{
          HashMap<String, int[]> index = StringIndexer.buildIndex(refSeq, counter.K);
          MappedRead mRead = counter.mapReadToRegion(seq, index, 0, refSeq.length());
          if(mRead != null){
            Alignment aln = aligner.align(seq,
                Arrays.copyOfRange(refSeq.getBytes(), mRead.start, mRead.end));
            scores.add(new RefScore(ref, (float)aln.score/fullScore, (float)aln.identity/aln.length,
                (float)(aln.end1 - aln.start1)/len));
          }
        }catch(Exception e){
          e.printStackTrace();
        }
      }
      if(scores.size() == 0){
        cluster.markedForDeletion = true;
        return;
      }
      Collections.sort(scores);
      RefScore max = scores.get(0);
      float maxScore = max.score;
      RefSubseq maxRef = max.ref;
      float maxIdentity = max.identity;
      cluster.refIDToScore = new TIntFloatHashMap();
      if(maxScore > 0){
        for(RefScore refScore: scores){
          if(maxScore - refScore.score <= delta){
            cluster.addRef(refScore.ref);
            cluster.refIDToScore.put(refScore.ref.refID, refScore.score);
          }
        }
        score = maxScore;
        identity = maxIdentity;
        if(changeCentroid){
          Reference reference = refAlignment.refAlns.get(maxRef.refID);
          cluster.centroidID = reference.name;
          cluster.centroidSeq = reference.seq.substring(maxRef.start, maxRef.end);
          cluster.centroidStart = maxRef.start;
          cluster.centroidEnd = maxRef.end;
        }
      }
    }

  }

  private void assignReferencesToClusters(LinkedList<Cluster>[] clusters){
    try{
      ArrayList<SWScore> SWTasks = new ArrayList<>();
      ExecutorService executorSW = Executors.newFixedThreadPool(threadNum);
      for(LinkedList<Cluster> clusters_i : clusters){
        for(Cluster cluster: clusters_i){
          if(cluster.contig == null || cluster.contig.seq.length == 0){
            cluster.markedForDeletion = true;
            continue;
          }
          SWScore task = new SWScore(cluster, true);
          SWTasks.add(task);
          executorSW.execute(task);
        }
      }
      executorSW.shutdown();
      executorSW.awaitTermination(10, TimeUnit.HOURS);
      for(int i = 0, j =0; i < clusters.length; i++){
        LinkedList<Cluster> clusters_i = clusters[i];
        for(Cluster cluster : clusters_i){
          if(cluster.markedForDeletion){
            continue;
          }
          SWScore task = SWTasks.get(j);
          j++;
          if(debug){
            if(cluster.refNum > 0){
              logger.printf("Cluster " + cluster.name + " read num = " + cluster.readNum +
                  " assigned to " + cluster.centroidID + " score = " + task.score +
                  " identity = " + task.identity + " coverage = " + task.coverage);
              for(TIntIterator it = cluster.refIDs.iterator(); it.hasNext(); ){
                String name = refAlignment.refAlns.get(it.next()).name;
                if(!name.equals(cluster.centroidID)){
                  logger.printf(" " + name);
                }
              }
              logger.println("");
            }else{
              logger.println("Can't find references for cluster " + cluster.name +
                  " with " + cluster.readNum + " reads");
            }
          }
        }
      }
      PairScoring[][] tasks = computePairCounts(clusters);
      for(PairScoring[] tasks_i: tasks){
        Arrays.sort(tasks_i);
        for(int j = 0; j < tasks_i.length; j++){
          Cluster cluster_ij = tasks_i[j].cluster;
          float[] countsNorm_ij = tasks_i[j].countsNorm;
          if(cluster_ij.markedForDeletion || cluster_ij.refNum == 0){
            continue;
          }
          TIntArrayList mergeList = new TIntArrayList();
          TObjectIntHashMap<String> centroidIDToCount = new TObjectIntHashMap<>();
          centroidIDToCount.put(cluster_ij.centroidID, cluster_ij.readNum);
          for(int k = j + 1; k < tasks_i.length; k++){
            Cluster cluster_ik = tasks_i[k].cluster;
            float[] countsNorm_ik = tasks_i[k].countsNorm;
            if(cluster_ik.markedForDeletion){// || cluster_ij.refNum != cluster_ik.refNum
              continue;
            }
            float sim = 0;
            for(int n = 0; n < countsNorm_ij.length; n++){
              sim += countsNorm_ij[n]*countsNorm_ik[n];
            }
            //if(cluster_ij.refIDs.containsAll(cluster_ik.refIDs)){
            if(sim >= similarityThreshold){
              mergeList.add(k);
              centroidIDToCount.adjustOrPutValue(cluster_ik.centroidID,
                  cluster_ik.readNum, cluster_ik.readNum);
            }
          }
          int maxCount = 0;
          String maxCentroidID = null;
          for(TObjectIntIterator<String> iter = centroidIDToCount.iterator(); iter.hasNext(); ){
            iter.advance();
            if(iter.value() > maxCount){
              maxCount = iter.value();
              maxCentroidID = iter.key();
            }
          }
          for(int k = 0; k < mergeList.size(); k++){
            Cluster cluster_ik = tasks_i[mergeList.getQuick(k)].cluster;
            cluster_ik.markedForDeletion = true;
            cluster_ij.readIDs.addAll(cluster_ik.readIDs);
            cluster_ij.readNum += cluster_ik.readNum;
            if(debug){
              logger.printf("Merging " + cluster_ik.name + " to " + cluster_ij.name);
            }
            if(!cluster_ij.centroidID.equals(maxCentroidID) &&
                cluster_ik.centroidID.equals(maxCentroidID)){
              cluster_ij.centroidID = maxCentroidID;
              cluster_ij.centroidSeq = cluster_ik.centroidSeq;
              cluster_ij.centroidStart = cluster_ik.centroidStart;
              cluster_ij.centroidEnd = cluster_ik.centroidEnd;
            }
            if(debug){
              logger.println(" centroid = " + cluster_ij.centroidID);
            }
          }
        }
      }
      LinkedList<Cluster> clustersToRebuildContigs = new LinkedList<>();
      ContigBuilder contigBuilder = new ContigBuilder(config, mappedData.mappedReads);
      for(LinkedList<Cluster> clusters_i : clusters){
        for(Iterator<Cluster> iter = clusters_i.iterator(); iter.hasNext();){
          Cluster cluster = iter.next();
          if(cluster.markedForDeletion){
            iter.remove();
          }else if(cluster.refNum > 0){
            clustersToRebuildContigs.add(cluster);
          }
        }
      }
      if(debug){
        logger.println("After merging: ");
        for(LinkedList<Cluster> clusters_i : clusters){
          for(Cluster cluster : clusters_i){
            logger.printf("Cluster " + cluster.name + " read num = " + cluster.readNum +
                " assigned to " + cluster.centroidID);
            for(TIntIterator it = cluster.refIDs.iterator(); it.hasNext(); ){
              String name = refAlignment.refAlns.get(it.next()).name;
              if(!name.equals(cluster.centroidID)){
                logger.printf(" " + name);
              }
            }
            logger.println("");
          }
        }
      }
      LinkedList<LinkedList<Cluster>> rebuiltContigs = contigBuilder.buildContigs(clustersToRebuildContigs);
      for(LinkedList<Cluster> clusters_i : clusters){
        int size = clusters_i.size();
        for(int j = 0; j < size; j++){
          Cluster cluster = clusters_i.pop();
          if(cluster.refNum > 0){
            clusters_i.addAll(rebuiltContigs.pop());
          }else{
            clusters_i.add(cluster);
          }
        }
      }
      executorSW = Executors.newFixedThreadPool(threadNum);
      for(LinkedList<Cluster> clusters_i : clusters){
        for(Cluster cluster: clusters_i){
          cluster.refNum = 0;
          cluster.refIDs.clear();
          cluster.refSubseqs.clear();
          SWScore task = new SWScore(cluster, false);
          executorSW.execute(task);
        }
      }
      executorSW.shutdown();
      executorSW.awaitTermination(10, TimeUnit.HOURS);
    }catch(Exception e){
      e.printStackTrace();
      System.exit(1);
    }
  }

}
