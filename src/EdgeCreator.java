/**
 * Created by Gennady on 16.09.2015.
 */
abstract class EdgeCreator{

  void createEdges(Cluster[] clusters){
    for(int i = 0; i < clusters.length; i++){
      Cluster cluster_i = clusters[i];
      for(int j = i + 1; j < clusters.length; j++){
        Cluster cluster_j = clusters[j];
        if(cluster_j.contig.alnStart < cluster_i.contig.alnEnd){
          createEdge(cluster_i, cluster_j);
        }else{
          break;
        }
      }
    }
  }

  void repairEdges(Cluster[] clusters){
    for(Cluster cluster: clusters){
      if(cluster.edgesNeedToBeRepaired){
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
        cluster.neighborsBackward.clear();
        cluster.neighborsForward.clear();
      }
    }
    for(int i = 0; i < clusters.length; i++){
      Cluster cluster_i = clusters[i];
      for(int j = i + 1; j < clusters.length; j++){
        Cluster cluster_j = clusters[j];
        if(!cluster_i.edgesNeedToBeRepaired && !cluster_j.edgesNeedToBeRepaired){
          continue;
        }
        if(cluster_j.contig.alnStart < cluster_i.contig.alnEnd){
          createEdge(cluster_i, cluster_j);
        }else{
          break;
        }
      }
    }
    for(Cluster cluster: clusters){
      cluster.edgesNeedToBeRepaired = false;
    }
  }

  abstract void createEdge(Cluster cluster_i, Cluster cluster_j);
}
