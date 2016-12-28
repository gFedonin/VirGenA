import gnu.trove.iterator.TIntIterator;
import org.jdom2.Document;

/**
 * Created by Gennady on 18.09.2015.
 */
class CommonReferenceAndReadEdgeCreator extends EdgeCreator{

  void createEdge(Cluster cluster_i, Cluster cluster_j){
    boolean haveCommonRef = false;
    if(cluster_i.refNum < cluster_j.refNum){
      for(TIntIterator iter = cluster_i.refIDs.iterator();iter.hasNext();){
        if(cluster_j.refIDs.contains(iter.next())){
          haveCommonRef = true;
          break;
        }
      }
    }else{
      for(TIntIterator iter = cluster_j.refIDs.iterator();iter.hasNext();){
        if(cluster_i.refIDs.contains(iter.next())){
          haveCommonRef = true;
          break;
        }
      }
    }
    float weight = (haveCommonRef)?1:0;
    int commonReadNum = 0;
    if(cluster_i.readNum < cluster_j.readNum){
      for(TIntIterator iter = cluster_i.readIDs.iterator(); iter.hasNext(); ){
        int id = iter.next();
        if(cluster_j.readIDs.contains(id)){
          commonReadNum++;
        }
      }
    }else{
      for(TIntIterator iter = cluster_j.readIDs.iterator(); iter.hasNext(); ){
        int id = iter.next();
        if(cluster_i.readIDs.contains(id)){
          commonReadNum++;
        }
      }
    }
    weight += (float)commonReadNum/(cluster_i.readNum + cluster_j.readNum - commonReadNum);
    weight /= 2;
    if(weight > 0){
      cluster_i.neighborsForward.add(cluster_j);
      cluster_j.neighborsBackward.add(cluster_i);
      cluster_i.weightsForward.add(weight);
      cluster_j.weightsBackward.add(weight);
    }
  }

}
