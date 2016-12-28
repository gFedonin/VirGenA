import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TIntArrayList;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;

/**
 * Created by Gennady on 17.04.2015.
 */
class Path extends Constants{

  int id;
  private LinkedList<Cluster> path;
  private TByteArrayList isNewVertex;
  float weight;
  int unvisitedNum;
  byte[] contig;
  byte[] contigAln;
  int[] coverage;
  int contigAlnStart;
  int contigAlnEnd;
  int refID;
  ArrayList<Contig>[] correctContigs;

  Path(){
    path = new LinkedList<>();
    isNewVertex = new TByteArrayList();
  }

  public Iterator<Cluster> iterator(){
    return path.iterator();
  }

  public int size(){
    return path.size();
  }

  void addVertex(Cluster vertex){
    if(path.size() > 0 && vertex.contig.alnStart > path.getFirst().contig.alnStart){
      System.out.println("Path not sorted");
      System.exit(1);
    }
    path.addFirst(vertex);
    vertex.weight = 0;
    if(vertex.visited){
      isNewVertex.add((byte)0);
    }else{
      unvisitedNum ++;
      isNewVertex.add((byte) 1);
      vertex.visited = true;
    }

  }

  void concatenateContigsRefID(){
    isNewVertex.reverse();
    TByteArrayList contigBuilder = new TByteArrayList();
    TIntArrayList cov = new TIntArrayList();
    TByteArrayList refIsCorrect = new TByteArrayList();
    Iterator<Cluster> iter = path.iterator();
    Cluster first = iter.next();
    contigAlnStart = first.contig.alnStart;
    contigAlnEnd = first.contig.alnEnd;
    contigBuilder.addAll(first.contig.seqAln);
    cov.addAll(first.contig.coverageAln);
    byte mark = first.refIDs.contains(refID)?(byte)1:0;
    for(int i = 0; i < first.contig.seqAln.length; i++){
      refIsCorrect.add(mark);
    }
    while(iter.hasNext()){
      Cluster curr = iter.next();
      mark = curr.refIDs.contains(refID)?(byte)1:0;
      if(curr.contig.alnEnd > contigAlnEnd){
        int n = curr.contig.alnStart - contigAlnStart;
        for(int m = 0; n < contigAlnEnd - contigAlnStart; n++, m++){
          if(mark == 1){
            if(refIsCorrect.getQuick(n) == 1){
              if(curr.contig.coverageAln[m] < cov.getQuick(n)){
                if(contigBuilder.getQuick(n) != curr.contig.seqAln[m]){
                  contigBuilder.setQuick(n, curr.contig.seqAln[m]);
                }
                cov.setQuick(n, curr.contig.coverageAln[m]);
              }
            }else{
              if(contigBuilder.getQuick(n) != curr.contig.seqAln[m]){
                contigBuilder.setQuick(n, curr.contig.seqAln[m]);
              }
              cov.setQuick(n, curr.contig.coverageAln[m]);
              refIsCorrect.setQuick(n, (byte)1);
            }
          }else{
            if(refIsCorrect.getQuick(n) == 0){
              if(curr.contig.coverageAln[m] < cov.getQuick(n)){
                if(contigBuilder.getQuick(n) != curr.contig.seqAln[m]){
                  contigBuilder.setQuick(n, curr.contig.seqAln[m]);
                }
                cov.setQuick(n, curr.contig.coverageAln[m]);
              }
            }
          }
        }
        for(int m = contigAlnEnd - curr.contig.alnStart;
            m < curr.contig.seqAln.length; m++){
          contigBuilder.add(curr.contig.seqAln[m]);
          cov.add(curr.contig.coverageAln[m]);
          refIsCorrect.add(mark);
        }
        contigAlnEnd = curr.contig.alnEnd;
      }else{
        int n = curr.contig.alnStart - contigAlnStart;
        for(int m = 0; n < curr.contig.alnEnd - contigAlnStart; n++, m++){
          if(mark == 1){
            if(refIsCorrect.getQuick(n) == 1){
              if(m >= curr.contig.coverageAln.length || n >= cov.size()){
                int a = 0;
              }
              if(curr.contig.coverageAln[m] < cov.getQuick(n)){
                if(contigBuilder.getQuick(n) != curr.contig.seqAln[m]){
                  contigBuilder.setQuick(n, curr.contig.seqAln[m]);
                }
                cov.setQuick(n, curr.contig.coverageAln[m]);
              }
            }else{
              if(contigBuilder.getQuick(n) != curr.contig.seqAln[m]){
                contigBuilder.setQuick(n, curr.contig.seqAln[m]);
              }
              cov.setQuick(n, curr.contig.coverageAln[m]);
              refIsCorrect.setQuick(n, (byte)1);
            }
          }else{
            if(refIsCorrect.getQuick(n) == 0){
              if(curr.contig.coverageAln[m] < cov.getQuick(n)){
                if(contigBuilder.getQuick(n) != curr.contig.seqAln[m]){
                  contigBuilder.setQuick(n, curr.contig.seqAln[m]);
                }
                cov.setQuick(n, curr.contig.coverageAln[m]);
              }
            }
          }
        }
      }
    }
    contigAln = contigBuilder.toArray();
    TByteArrayList finContig = new TByteArrayList();
    TIntArrayList finCov = new TIntArrayList();
    for(int i = 0; i < contigBuilder.size(); i++){
      if(contigBuilder.getQuick(i) != GAP){
        finContig.add(contigBuilder.getQuick(i));
        finCov.add(cov.getQuick(i));
      }
    }
    contig = finContig.toArray();
    coverage = finCov.toArray();
  }

}