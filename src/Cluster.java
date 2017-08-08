import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TIntFloatHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.TIntHashSet;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by Gennady on 12.09.2015.
 */
class Cluster extends Constants implements Comparable{

  int id;
  int windowID;
  int start;
  int end;
  String name;
  Contig contig;
  int centroidStart;
  int centroidEnd;
  String centroidID;
  String centroidSeq;
  TIntHashSet refIDs;
  TIntFloatHashMap refIDToScore;
  ArrayList<RefSubseq> refSubseqs;
  TIntObjectHashMap<MappedRead> mReads;
  int readNum;
  int refNum;
  TIntHashSet readIDs;
  TIntHashSet readsToRemap;
  boolean markedForDeletion;
  boolean rebuildContig;
  HashMap<String, int[]> centroidIndex;
  ArrayList<Cluster> neighborsForward;
  ArrayList<Cluster> neighborsBackward;
  TFloatArrayList weightsForward;
  TFloatArrayList weightsBackward;
  float[] vertexVals;
  boolean visited;
  boolean startsPath;
  boolean endsPath;
  float weight;
  boolean edgesNeedToBeRepaired;

  Cluster(int id, int start, int end, String centroidID){
    this.id = id;
    this.start = start;
    this.end = end;
    this.centroidID = centroidID;
    refIDs = new TIntHashSet();
    refSubseqs = new ArrayList<>();
    readIDs = new TIntHashSet();
    readsToRemap = new TIntHashSet();
    name = String.valueOf(start) + '_' + end + '_' + id;
    neighborsForward = new ArrayList<>();
    neighborsBackward = new ArrayList<>();
    weightsForward = new TFloatArrayList();
    weightsBackward = new TFloatArrayList();
  }

  Cluster(Cluster cluster, byte[] contigSeq,
                    int[] coverage){
    id = cluster.id;
    windowID = cluster.windowID;
    start = cluster.start;
    end = cluster.end;
    centroidID = cluster.centroidID;
    refIDs = cluster.refIDs;
    refSubseqs = cluster.refSubseqs;
    readIDs = new TIntHashSet();
    readsToRemap = new TIntHashSet();
    //readIDs = new MyBitSet();
    name = String.valueOf(start) + '_' + String.valueOf(end) + '_' + String.valueOf(id);
    centroidSeq = cluster.centroidSeq;
    centroidStart = cluster.centroidStart;
    centroidEnd = cluster.centroidEnd;
    contig = new Contig(cluster.contig);
    contig.seq = contigSeq;
    contig.coverage = coverage;
    mReads = new TIntObjectHashMap<>();
    for(TIntObjectIterator<MappedRead> it = cluster.mReads.iterator();it.hasNext();){
      it.advance();
      int id = it.key();
      MappedRead mRead = it.value();
      if(mRead.start + mRead.aln.start2 >= contig.centroidStart &&
          mRead.start + mRead.aln.end2 <= contig.centroidEnd){
        readIDs.add(id);
        //readIDs.set(id);
        mReads.put(id, mRead);
      }
    }
    readNum = readIDs.size();
    refNum = refIDs.size();
    neighborsForward = new ArrayList<>();
    neighborsBackward = new ArrayList<>();
    weightsForward = new TFloatArrayList();
    weightsBackward = new TFloatArrayList();
  }


  void buildCentroidIndex(int K){
    centroidIndex = StringIndexer.buildIndex(centroidSeq, K);
  }

  @Override
  public int compareTo(Object o){
    Cluster cluster = (Cluster) o;
    if(contig.alnStart == cluster.contig.alnStart){
      return cluster.contig.alnEnd - contig.alnEnd;
    }
    return contig.alnStart - cluster.contig.alnStart;
  }

  void addRef(RefSubseq ref){
    refNum ++;
    refIDs.add(ref.refID);
    refSubseqs.add(ref);
  }

  void addRead(int readID){
    readNum ++;
    readIDs.add(readID);
  }

}
