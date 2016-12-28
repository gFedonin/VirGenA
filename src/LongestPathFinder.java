import gnu.trove.iterator.TFloatIterator;
import gnu.trove.iterator.TIntIterator;
import org.jdom2.Document;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;

/**
 * Created by Gennady on 17.09.2015.
 */
class LongestPathFinder{

  private ReferenceAlignment refAlignment;
  Logger logger;
  private boolean debug;

  LongestPathFinder(Document document){
    refAlignment = ReferenceAlignment.getInstance(document);
    logger = Logger.getInstance(document);
    debug = Boolean.parseBoolean(
        document.getRootElement().getChild("ReferenceSelector").getChildText("Debug"));
  }

  ArrayList<Path> findLongestPaths(Cluster[] vertices){

    ArrayList<Path> longestPaths = new ArrayList<>();
    LinkedList<Cluster> stack = new LinkedList<>();
    for(Cluster cluster : vertices){
      cluster.vertexVals = new float[refAlignment.refAlns.size()];
      if(!cluster.visited){
        topologicalSort(cluster, stack);
      }
    }

    int clusterNum = stack.size();
    int visitedNum = 0;
    for(Cluster vertex : stack){
      vertex.visited = false;
      if(vertex.neighborsBackward.size() == 0){
        vertex.startsPath = true;
      }
      if(vertex.neighborsForward.size() == 0){
        vertex.endsPath = true;
      }
    }

    int i = 0;
    while(visitedNum != clusterNum){
      for(Cluster vertex : stack){
        for(int j = 0; j < vertex.vertexVals.length; j++){
          vertex.vertexVals[j] = 0;
        }
      }
      for(Cluster vertex : stack){
        Iterator<Cluster> iterator = vertex.neighborsForward.iterator();
        TFloatIterator weightIter = vertex.weightsForward.iterator();
        while(iterator.hasNext()){
          Cluster neighbor = iterator.next();
          float weight = weightIter.next() + vertex.weight;
          for(TIntIterator iter = vertex.refIDs.iterator(); iter.hasNext();){
            int refID = iter.next();
            if(neighbor.refIDs.contains(refID)){
              if(weight + vertex.vertexVals[refID] > neighbor.vertexVals[refID]){
                neighbor.vertexVals[refID] = weight + vertex.vertexVals[refID];
              }
            }
          }
        }
      }
      for(Cluster vertex: vertices){
        for(TIntIterator iter = vertex.refIDs.iterator();iter.hasNext();){
          int refID = iter.next();
          boolean endPath = true;
          for(Cluster neighbor: vertex.neighborsForward){
            if(neighbor.refIDs.contains(refID)){
              endPath = false;
              break;
            }
          }
          if(endPath){
            vertex.vertexVals[refID] += vertex.weight;
          }
        }
      }
      Path maxPath = backwardTraverse(vertices);
      if(maxPath.unvisitedNum == 0){
        break;
      }
      visitedNum += maxPath.unvisitedNum;
      maxPath.id = i;
      i++;
      //maxPath.concatenateContigsVertexMark();
      maxPath.concatenateContigsRefID();
      longestPaths.add(maxPath);
    }

    return longestPaths;
  }

  private Path backwardTraverse(Cluster[] vertices){
    Cluster currVertex = null;
    float maxVal = Float.NEGATIVE_INFINITY;
    int maxRef = 0;
    for(Cluster cluster : vertices){
      float maxRefVal = 0;
      int maxRefId = 0;
      for(int j = 0; j < cluster.vertexVals.length; j++){
        if(cluster.vertexVals[j] > maxRefVal){
          maxRefVal = cluster.vertexVals[j];
          maxRefId = j;
        }
      }
      if(maxRefVal > maxVal){
        currVertex = cluster;
        maxRef = maxRefId;
        maxVal = maxRefVal;
      }
    }
    Path path = new Path();
    path.weight = maxVal;
    path.refID = maxRef;
    path.addVertex(currVertex);
    while(currVertex.neighborsBackward.size() > 0){
      TFloatIterator iterW = currVertex.weightsBackward.iterator();
      Cluster backNeighbor = null;
      boolean visited = false;
      float max = Float.NEGATIVE_INFINITY;
      for(Cluster vertex: currVertex.neighborsBackward){
        float weight = iterW.next() + vertex.weight;
        if(!vertex.visited){
          if(vertex.refIDs.contains(maxRef)){
            if(vertex.vertexVals[maxRef] + weight > max){
              max = vertex.vertexVals[maxRef] + weight;
              backNeighbor = vertex;
              visited = false;
            }else{
              if(vertex.vertexVals[maxRef] + weight == max){
                if(visited){
                  backNeighbor = vertex;
                  visited = false;
                }
              }
            }
          }
        }else{
          if(vertex.refIDs.contains(maxRef)){
            if(vertex.vertexVals[maxRef] + weight > max){
              max = vertex.vertexVals[maxRef] + weight;
              backNeighbor = vertex;
              visited = true;
            }
          }
        }
      }
      if(backNeighbor == null){
        break;
      }
      currVertex = backNeighbor;
      path.addVertex(backNeighbor);
    }
    if(debug){
      logger.printf("Found path for ref %s of weight %1.2f\n",
          refAlignment.refAlns.get(maxRef).name, path.weight);
    }
    return path;
  }

  private void topologicalSort(Cluster vertex, LinkedList<Cluster> stack){
    // Mark the current node as visited
    vertex.visited = true;
    // Recur for all the vertices adjacent to this vertex
    for(Cluster neighbor : vertex.neighborsForward){
      if(!neighbor.visited){
        topologicalSort(neighbor, stack);
      }
    }
    // Push current vertex to stack which stores topological sort
    stack.push(vertex);
  }

}
