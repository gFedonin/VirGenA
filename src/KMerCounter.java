import org.jdom2.Document;
import org.jdom2.Element;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by Gennady on 06.12.2015.
 */
class KMerCounter extends KMerCounterBase{

  private static KMerCounter instance;

  private KMerCounter(Document document){
    try{
      Element element = document.getRootElement().getChild("Mapper");
      DataReader dataReader = DataReader.getInstance(document);
      pValue = Float.parseFloat(element.getChildText("pValue"));
//      randomReadsCountsPath = element.getChildText("RandomModelPath");
//      if(randomReadsCountsPath.isEmpty()){
        K = Integer.parseInt(element.getChildText("K"));
        element = element.getChild("RandomModelParameters");
//        String outPath = document.getRootElement().getChildText("OutPath");
        Logger logger = Logger.getInstance(document);
//        logger.println("No path to random read model found in the config file for <Mapper>. Creating new model with the " +
//            "parameters given in the config file in <RandomModelParameters> and K = " + scoreK);
        logger.println("Creating random reads model from reference with the parameters given in the config file in <RandomModelParameters>");
        RandomModel model = new RandomModel(document);
//        int minReadLen = Integer.parseInt(element.getChildText("MinReadLen"));
//        int maxReadLen = Integer.parseInt(element.getChildText("MaxReadLen"));
        int step = Integer.parseInt(element.getChildText("Step"));
        int readNum = Integer.parseInt(element.getChildText("ReadNum"));
//        randomReadsCountsPath = outPath + "random_model.rm";
//        model.printModel(readNum, minReadLen, maxReadLen, step, randomReadsCountsPath);
//        logger.println("Random model is created and saved in " + randomReadsCountsPath);
//        logger.println("You can reuse it in future by providing the path to random_model.rm in <RandomModelPath> in " +
//            "config file in <Mapper> section. In case then such path is provided, the model parameters in " +
//            "<RandomModelParameters> section and the K in <Mapper> section are ignored.");
      model.genModel(readNum, dataReader.minReadLen, dataReader.maxReadLen, step);
//      }
//      readRandomModel();
      computeCuttofs(dataReader.minReadLen, dataReader.maxReadLen, step, model.counts);
    }catch(Exception e){
      e.printStackTrace();
    }
  }

  KMerCounter(){}

  private KMerCounter(int K, float coef){
    this.K = K;
    highCoef = coef;
    lowCoef = 1.0f/coef;
  }

  static KMerCounter getInstance(Document document){
    if(instance == null){
      instance = new KMerCounter(document);
    }
    return instance;
  }

  static KMerCounter getInstance(int K, float coef){
    if(instance == null){
      instance = new KMerCounter(K, coef);
    }
    return instance;
  }


  private Location getInitState(int[] contigEnds, int readLen,
                                LongHit[] longHits, boolean[] isAdjacent){
    Location res = new Location();
    int totalHit = longHits.length;
    LongHit h_i = null;
    outer:
    for(int i = 0; i < totalHit; i++){
      h_i = longHits[i];
//    int startMax = 0;
//    int endMax = 0;
      for(int j = 0, start = 0; j < contigEnds.length; j++){
        int end = contigEnds[j];
        if(start <= h_i.genomePos && h_i.genomePos < end){
          res.start = Math.max(start, h_i.genomePos - h_i.readPos*3/2);
          res.end = Math.min(h_i.genomePos + h_i.len + K + (readLen - h_i.readPos - h_i.len - K)*3/2, end);
          if(h_i.genomePos + h_i.len + K <= res.end){
            break outer;
          }
        }
//      if(locEnd - locStart > endMax - startMax){
//        endMax = locEnd;
//        startMax = locStart;
//      }
        start = end;
      }
    }
//    res.start = startMax;
//    res.end = endMax;
    res.startIndex = 0;
    int count = h_i.len;
    int j;
    for(j = 1; j < totalHit ; j++){
      LongHit h_j = longHits[j];
      if(h_j.genomePos + K + h_j.len > res.end){
        break;
      }
      count += h_j.len;
      if(isAdjacent[j - 1]){
        count ++;
      }
    }
    res.endIndex = j - 1;
    res.count = count;
    return res;
  }


  MappedRead mapReadToRegion(byte[] read, HashMap<String, int[]> index, int from, int to){
    int readLen = read.length;
    String readStr = new String(read);
    LongHit[] longHits = getHits(readStr, index);
    int totalHit = longHits.length;
    if(totalHit == 0){
      return null;
    }
    boolean[] isAdjacent = concordanceArray(longHits);
    Location state = getInitState(from, to, readLen, longHits, isAdjacent);
    MappedRead res = new MappedRead(state.start, state.end, state.count);
    for(int i = 1; i < totalHit; i++){
      computeCount(from, to, readLen, longHits, i, state, isAdjacent);
      if(state.count > res.count){
        res.count = state.count;
        res.start = state.start;
        res.end = state.end;
      }
    }
    res.seq = read;
    return res;
  }

  private LongHit[] splitLongHits(LongHit[] longHits, int[] contigEnds){
    int contigEnd = contigEnds[0];
    LongHit hit = longHits[0];
    ArrayList<LongHit> splitted = new ArrayList<>();
    int i = 0, j = 0;
    while(true){
      if(hit.genomePos < contigEnd){
        if(hit.genomePos + K + hit.len > contigEnd){
          // hit intersects with the next contig
          if(contigEnd - hit.genomePos >= K){
            //first part is long enough -> create new hit
            LongHit newHit = new LongHit(hit.genomePos, hit.readPos);
            newHit.len = contigEnd - hit.genomePos - K;
            splitted.add(newHit);
          }
          if(hit.genomePos + K + hit.len - contigEnd >= K){
            //second part is long enough -> truncate hit from left to the end of contig
            hit.readPos += contigEnd - hit.genomePos;
            hit.len -= contigEnd - hit.genomePos;
            hit.genomePos = contigEnd;
          }else{
            j ++;
            if(j < longHits.length){
              hit = longHits[j];
            }else{
              break;
            }
          }
        }else{
          splitted.add(hit);
          j ++;
          if(j < longHits.length){
            hit = longHits[j];
          }else{
            break;
          }
        }
      }else{
        i ++;
        contigEnd = contigEnds[i];
      }
    }
    return splitted.toArray(new LongHit[splitted.size()]);
  }

  MappedRead[] moveWindow(LongHit[] longHits, boolean[] isAdjacent, int[] contigEnds, Location location, int readLen){
  int totalHit = longHits.length;
    int cutoff = cutoffs[readLen];
    ArrayList<MappedRead> bestPositions = new ArrayList<>();
    if(location.count > cutoff){
      bestPositions.add(new MappedRead(location.start, location.end, location.count));
    }
    for(int i = 1; i < totalHit; i++){
//      int count = matches[location.startIndex].len;
//      for(int j = location.startIndex + 1; j <= location.endIndex; j++){
//        count += matches[j].len;
//        if(isAdjacent[j - 1]){
//          count ++;
//        }
//      }
//      if(count != location.count){
//        System.out.println("!!!");
//        int a = 0;
//      }
//      if(i == 151){
//        int a = 0;
//      }
      computeCount(contigEnds, readLen, longHits, i, location, isAdjacent);
//      if(location.count < 0){
//        int a = 0;
//      }
      if(location.count > cutoff){
        bestPositions.add(new MappedRead(location.start, location.end, location.count));
      }
    }
    if(bestPositions.size() == 0){
      return null;
    }
    return pruneWindows(bestPositions);
  }

  MappedRead[] getNBestRegions(byte[] read, HashMap<String, int[]> index, int[] contigEnds){
    String readStr = new String(read);
    LongHit[] longHits = getHits(readStr, index);
    if(longHits.length == 0){
      return null;
    }
    longHits = splitLongHits(longHits, contigEnds);
    int totalHit = longHits.length;
    if(totalHit == 0){
      return null;
    }
    boolean[] isAdjacent = concordanceArray(longHits);
    Location location = getInitState(contigEnds, read.length, longHits, isAdjacent);
    return moveWindow(longHits, isAdjacent, contigEnds, location, read.length);
  }

}
