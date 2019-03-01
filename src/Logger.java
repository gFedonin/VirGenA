import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.TIntHashSet;
import org.jdom2.Document;
import org.jdom2.Element;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * Date: 09.08.15
 */
public class Logger extends Constants{

  String outPath;
  private int minContigLen;
  private PrintStream writer;

  private static Logger instance;

  private Logger(Document document){
    Element referenceSelectorElement = document.getRootElement().getChild("ReferenceSelector");
    outPath = document.getRootElement().getChildText("OutPath");
    minContigLen = Integer.parseInt(referenceSelectorElement.getChildText("MinContigLength"));
    String logPath = outPath + "log.txt";
    try{
      File outFolder = new File(outPath);
      if(!outFolder.exists()){
        outFolder.mkdirs();
      }
      writer = new PrintStream(new FileOutputStream(logPath), true);
    }catch(Exception e){
      e.printStackTrace();
    }
  }

  static Logger getInstance(Document document){
    String outPath = document.getRootElement().getChildText("OutPath");
    if(instance == null){
      instance = new Logger(document);
    }else if(!instance.outPath.equals(outPath)){
      try{
        instance.writer.close();
      }catch(Exception e){
        e.printStackTrace();
      }
      instance = new Logger(document);
    }
    return instance;
  }

  void printf(String s, Object... args){
    try{
      writer.printf(s, args);
    }catch(Exception e){
      e.printStackTrace();
    }
  }

  void println(String s){
    try{
      writer.println(s);
    }catch(Exception e){
      e.printStackTrace();
    }
  }

  void close(){
    try{
      writer.close();
    }catch(Exception e){
      e.printStackTrace();
    }
  }

  void printTotalReadNumInGroups(TIntArrayList[] readGroups){
    TIntHashSet readIDs = new TIntHashSet();
    for(TIntArrayList readGroup: readGroups){
      readIDs.addAll(readGroup);
    }
    printf("Total reads in all groups: %d\n", readIDs.size());
  }

  void printTotalReadNumInClusters(LinkedList<Cluster>[] clusters){
    TIntHashSet readIDs = new TIntHashSet();
    int clusterNum = 0;
    for(LinkedList<Cluster> list: clusters){
      for(Cluster cluster: list){
        clusterNum ++;
        readIDs.addAll(cluster.readIDs);
      }
    }
    printf("Total reads in all %d vertices: %d\n", clusterNum, readIDs.size());
  }

  void printTotalReadNumInClusters(List<Cluster> clusters){
    TIntHashSet readIDs = new TIntHashSet();
    int clusterNum = 0;
    for(Cluster cluster: clusters){
      clusterNum ++;
      readIDs.addAll(cluster.readIDs);
    }
    printf("Total reads in all %d vertices: %d\n", clusterNum, readIDs.size());
  }


  void printTotalReadNumInGraph(Cluster[] clusters){
    TIntHashSet readIDs = new TIntHashSet();
    for(Cluster cluster: clusters){
      readIDs.addAll(cluster.readIDs);
    }
    printf("Total reads in graph: %d\n", readIDs.size());
  }


  void printTotalReadNumInPaths(ArrayList<Path> longestPaths){
    TIntHashSet readIDs = new TIntHashSet();
    for(Path path: longestPaths){
      for(Iterator<Cluster> iter = path.iterator(); iter.hasNext();){
        Cluster vertex = iter.next();
        readIDs.addAll(vertex.readIDs);
      }
    }
    printf("Total reads in graph paths: %d\n", readIDs.size());
  }

  void printTotalReadNumInLongPaths(ArrayList<Path> longestPaths){
    TIntHashSet readIDs = new TIntHashSet();
    for(Path path: longestPaths){
      if(path.contig.length >= minContigLen){
        for(Iterator<Cluster> iter = path.iterator(); iter.hasNext();){
          Cluster vertex = iter.next();
          readIDs.addAll(vertex.readIDs);
        }
      }
    }
    printf("Total reads in graph paths: %d\n", readIDs.size());
  }

  void printClusterStats(LinkedList<Cluster>[] clusters, int windowLen, int winStep){
    println("HSP\tCluster num\tCluster read num\tContig lengths\tAv read identities");
    for(int i = 0; i < clusters.length; i++){
      LinkedList<Cluster> clusterList = clusters[i];
      int start = i*winStep;
      int end = start + windowLen;
      StringBuilder builder = new StringBuilder();
      int clusterNum = 0;
      StringBuilder sizes = new StringBuilder();
      StringBuilder contigs = new StringBuilder();
      StringBuilder identities= new StringBuilder();
      for(Cluster cluster: clusterList){
        if(cluster.readNum == 0){
          continue;
        }
        clusterNum ++;
        sizes.append(cluster.readNum).append(',');
        if(cluster.contig == null){
          println("Contig == null");
          System.exit(1);
        }
        if(cluster.contig.seq == null){
          println("Contig.seqB = null");
          println("rebuildContig = " + cluster.rebuildContig);
          println("readNum = " + cluster.readNum);
          println("refNum = " + cluster.refNum);
          println("refIDs.size = " + cluster.refIDs.size());
          println("contig.alnStart = " + cluster.contig.alnStart);
          println("contig.alnEnd = " + cluster.contig.alnEnd);
          println("centroidID = " + cluster.centroidID);
          System.exit(1);
        }
        contigs.append(cluster.contig.seq.length).append(',');
        float avIdentity = 0;
        for(MappedRead mRead : cluster.mReads.valueCollection()){
          avIdentity += (float) mRead.aln.identity/mRead.aln.length;
        }
        identities.append(avIdentity/cluster.readNum).append(',');
      }
      builder.append(start).append('_').append(end).append('\t');
      builder.append(clusterNum).append('\t');
      builder.append(sizes.toString()).append('\t');
      builder.append(contigs.toString()).append('\t');
      builder.append(identities.toString());
      println(builder.toString());
    }
  }

  void printClusterStats1(LinkedList<Cluster>[] clusters, int windowLen, int winStep, int refID1, int refID2){
    println("HSP\tCluster num\tReads with fragment\tReads with ref1\tReads with ref2\tCluster read num");
    for(int i = 0; i < clusters.length; i++){
      LinkedList<Cluster> clusterList = clusters[i];
      int start = i*winStep;
      int end = start + windowLen;
      StringBuilder builder = new StringBuilder();
      int clusterNum = 0;
      int totalReads = 0;
      int readsWithReferences = 0;
      int readsWithRef1 = 0;
      int readsWithRef2 = 0;
      StringBuilder sizes = new StringBuilder();
      for(Cluster cluster: clusterList){
        if(cluster.readNum == 0){
          continue;
        }
        clusterNum ++;
        totalReads += cluster.readNum;
        if(cluster.refNum > 0){
          readsWithReferences += cluster.readNum;
          if(cluster.refIDs.contains(refID1)){
            readsWithRef1 += cluster.readNum;
            if(cluster.refIDs.contains(refID2)){
              readsWithRef2 += cluster.readNum;
              sizes.append(cluster.readNum).append("(12),");
            }else{
              sizes.append(cluster.readNum).append("(1),");
            }
          }else if(cluster.refIDs.contains(refID2)){
            readsWithRef2 += cluster.readNum;
            sizes.append(cluster.readNum).append("(2),");
          }else{
            sizes.append(cluster.readNum).append("(3),");
          }
        }else{
          sizes.append(cluster.readNum).append("(-)").append(',');
        }
      }
      builder.append(start).append('_').append(end).append('\t');
      builder.append(clusterNum).append('\t').append((float)readsWithReferences/totalReads);
      builder.append('\t').append((float)readsWithRef1/totalReads);
      builder.append('\t').append((float)readsWithRef2/totalReads);
      builder.append('\t').append(sizes.toString());
      println(builder.toString());
    }
  }

  void printReadNumWithNoRef(List<Cluster> clusters, int windowLen, int winStep, int windowNum){
    int[] readNumNoRef = new int[windowNum];
    int[] readNumTotal = new int[windowNum];
    for(Cluster cluster: clusters){
      if(cluster.refNum == 0){
        readNumNoRef[cluster.windowID] += cluster.readNum;
      }
      readNumTotal[cluster.windowID] += cluster.readNum;
    }
    println("Read num in clusters with no fragment:");
    for(int i = 0; i < windowNum; i++){
      int start = i*winStep;
      int end = start + windowLen;
      printf("%d_%d %d (%1.2f)\n", start, end, readNumNoRef[i], (float) readNumNoRef[i]/readNumTotal[i]);
    }
  }

  void printClusterStats(ArrayList<Cluster> clustersArr, int windowLen, int winStep, int windowNum){
    TIntObjectHashMap<LinkedList<Cluster>> clustersMap = new TIntObjectHashMap<>();
    for(Cluster cluster: clustersArr){
      LinkedList<Cluster> list = clustersMap.get(cluster.windowID);
      if(list == null){
        list = new LinkedList<>();
        clustersMap.put(cluster.windowID, list);
      }
      list.add(cluster);
    }
    LinkedList<Cluster>[] clusters = new LinkedList[windowNum];
    for(TIntObjectIterator<LinkedList<Cluster>> iter = clustersMap.iterator(); iter.hasNext();){
      iter.advance();
      clusters[iter.key()] = iter.value();
    }
    println("HSP\tCluster num\tCluster read num\tContig lengths\tAv read identities");
    for(int i = 0; i < windowNum; i++){
      LinkedList<Cluster> clusterList = clusters[i];
      if(clusterList == null){
        continue;
      }
      int start = i*winStep;
      int end = start + windowLen;
      StringBuilder builder = new StringBuilder();
      int clusterNum = 0;
      StringBuilder sizes = new StringBuilder();
      StringBuilder contigs = new StringBuilder();
      StringBuilder identities= new StringBuilder();
      for(Cluster cluster: clusterList){
        clusterNum ++;
        sizes.append(cluster.readNum).append(',');
        if(cluster.contig == null){
          contigs.append('-').append(',');
          identities.append('-').append(',');
          continue;
        }
        if(cluster.contig.seq == null){
          println("Contig.seqB = null");
          println("rebuildContig = " + cluster.rebuildContig);
          println("readNum = " + cluster.readNum);
          println("refNum = " + cluster.refNum);
          println("refIDs.size = " + cluster.refIDs.size());
          println("contig.alnStart = " + cluster.contig.alnStart);
          println("contig.alnEnd = " + cluster.contig.alnEnd);
          println("centroidID = " + cluster.centroidID);
          System.exit(1);
        }
        contigs.append(cluster.contig.seq.length).append(',');
        float avIdentity = 0;
        for(MappedRead mRead : cluster.mReads.valueCollection()){
          avIdentity += (float) mRead.aln.identity/mRead.aln.length;
        }
        identities.append(avIdentity/cluster.readNum).append(',');
      }
      builder.append(start).append('_').append(end).append('\t');
      builder.append(clusterNum).append('\t');
      builder.append(sizes.toString()).append('\t');
      builder.append(contigs.toString()).append('\t');
      builder.append(identities.toString());
      println(builder.toString());
    }
  }

  void printVertexStats(Cluster[] vertices, int windowLen, int winStep, int windowNum){
    TIntObjectHashMap<LinkedList<Cluster>> clustersMap = new TIntObjectHashMap<>();
    for(Cluster cluster: vertices){
      LinkedList<Cluster> list = clustersMap.get(cluster.windowID);
      if(list == null){
        list = new LinkedList<>();
        clustersMap.put(cluster.windowID, list);
      }
      list.add(cluster);
    }
    LinkedList<Cluster>[] clusters = new LinkedList[windowNum];
    for(TIntObjectIterator<LinkedList<Cluster>> iter = clustersMap.iterator(); iter.hasNext();){
      iter.advance();
      clusters[iter.key()] = iter.value();
    }
    println("Cluster (contig) num\tCluster read num\tContig lengths\tAv read identities");
    for(int i = 0; i < windowNum; i++){
      LinkedList<Cluster> clusterList = clusters[i];
      if(clusterList == null){
        continue;
      }
      int start = i*winStep;
      int end = start + windowLen;
      StringBuilder builder = new StringBuilder();
      int clusterNum = 0;
      StringBuilder sizes = new StringBuilder();
      StringBuilder contigs = new StringBuilder();
      StringBuilder identities= new StringBuilder();
      for(Cluster cluster: clusterList){
        if(cluster.readNum > 0){
          clusterNum ++;
          sizes.append(cluster.readNum).append(',');
          if(cluster.contig == null){
            println("Contig == null");
            System.exit(1);
          }
          if(cluster.contig.seq == null){
            println("Contig.seqB = null");
            println("rebuildContig = " + cluster.rebuildContig);
            println("readNum = " + cluster.readNum);
            println("refNum = " + cluster.refNum);
            println("refIDs.size = " + cluster.refIDs.size());
            println("contig.alnStart = " + cluster.contig.alnStart);
            println("contig.alnEnd = " + cluster.contig.alnEnd);
            println("centroidID = " + cluster.centroidID);
            System.exit(1);
          }
          contigs.append(cluster.contig.seq.length).append(',');
          float avIdentity = 0;
          for(MappedRead mRead : cluster.mReads.valueCollection()){
            avIdentity += (float) mRead.aln.identity/mRead.aln.length;
          }
          identities.append(avIdentity/cluster.readNum).append(',');
        }
      }
      builder.append(start).append('_').append(end).append('\t');
      builder.append(clusterNum).append('\t');
      builder.append(sizes.toString()).append('\t');
      builder.append(contigs.toString()).append('\t');
      builder.append(identities.toString());
      println(builder.toString());
    }
  }

}
