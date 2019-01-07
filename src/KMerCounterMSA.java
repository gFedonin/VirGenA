import org.jdom2.Document;
import org.jdom2.Element;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * Date: 03.01.15
 */
class KMerCounterMSA extends KMerCounterBase{

  private static KMerCounterMSA instance;

  private KMerCounterMSA(Document document){
    Element element = document.getRootElement().getChild("ReferenceSelector").getChild("MapperToMSA");
    randomReadsCountsPath = element.getChildText("RandomModelPath");
    pValue = Float.parseFloat(element.getChildText("pValue"));
    try {
      cutoffs = readRandomModel();
    }catch (Exception e){
      e.printStackTrace();
    }
  }

  private KMerCounterMSA(int K, float coef){
    this.K = K;
    this.highCoef = coef;
    this.lowCoef = 1/coef;
  }

  static KMerCounterMSA getInstance(Document document){
    if(instance == null){
      instance = new KMerCounterMSA(document);
    }
    return instance;
  }

  static KMerCounterMSA getInstance(int K, float coef){
    if(instance == null){
      instance = new KMerCounterMSA(K, coef);
    }
    return instance;
  }

  private LongHit[] mergeHits(AlnPos[][] hits){
    ArrayList<LongHit> longHits = new ArrayList<>();
    ArrayList<LongHit> currSet = new ArrayList<>();
    ArrayList<LongHit> nextSet = new ArrayList<>();
    int firstMatch = 0;
    do{
      if(firstMatch == hits.length){
        return new LongHit[0];
      }
      for(int j = 0; j < hits[firstMatch].length; j++){
        AlnPos pos_j = hits[firstMatch][j];
        LongHit longHit = new LongHit(pos_j.pos, firstMatch);
        currSet.add(longHit);
        longHits.add(longHit);
      }
      firstMatch ++;
    }while(currSet.isEmpty());
    for(int i = firstMatch; i < hits.length; i++){
      AlnPos[] hits_i = hits[i];
      if(hits_i.length == 0){
        currSet = new ArrayList<>();
        nextSet = new ArrayList<>();
        continue;
      }
      if(currSet.size() == 0){
        for(AlnPos pos_j : hits_i){
          LongHit longHit = new LongHit(pos_j.pos, i);
          currSet.add(longHit);
          longHits.add(longHit);
        }
        continue;
      }
      Iterator<LongHit> iteratorSet = currSet.iterator();
      LongHit lHit = iteratorSet.next();
      int j = 0;
      AlnPos rHit = hits_i[j];
      while(true){
        // no intersection
        if(lHit.genomePos + lHit.len + 1 < rHit.pos){
          //longHits.add(lHit);
          if(iteratorSet.hasNext()){
            lHit = iteratorSet.next();
          }else{
            LongHit longHit = new LongHit(rHit.pos, i);
            nextSet.add(longHit);
            longHits.add(longHit);
            j++;
            break;
          }
        }else if(lHit.genomePos + lHit.len + 1 > rHit.pos){
          j++;
          LongHit longHit = new LongHit(rHit.pos, i);
          nextSet.add(longHit);
          longHits.add(longHit);
          if(j < hits_i.length){
            rHit = hits_i[j];
          }else{
            break;
          }
        }else{
          nextSet.add(lHit);
          lHit.len++;
          j++;
          if(iteratorSet.hasNext()){
            lHit = iteratorSet.next();
          }else{
            break;
          }
          if(j < hits_i.length){
            rHit = hits_i[j];
          }else{
            break;
          }
        }
      }
      while(j < hits_i.length){
        rHit = hits_i[j];
        LongHit longHit = new LongHit(rHit.pos, i);
        nextSet.add(longHit);
        longHits.add(longHit);
        j ++;
      }
      currSet = nextSet;
      nextSet = new ArrayList<>();
    }
    TreeSet<LongHit> sortedSet = new TreeSet<>(longHits);
    ArrayList<LongHit> res = new ArrayList<>();
    LongHit prev = sortedSet.pollFirst();
    while(!sortedSet.isEmpty()){
      LongHit curr = sortedSet.pollFirst();
      if(prev.genomePos + prev.len + K - 1 >= curr.genomePos){
        if(prev.len < curr.len){
          prev.len = curr.genomePos - prev.genomePos - K;
          if(prev.len >= 0){
            res.add(prev);
          }
          prev = curr;
        }else{
          int lenDiff = prev.genomePos + prev.len + K - curr.genomePos;
          curr.len -= lenDiff;
          if(curr.len >= 0){
            curr.genomePos += lenDiff;
            curr.readPos += lenDiff;
            sortedSet.add(curr);
          }
        }
      }else{
        res.add(prev);
        prev = curr;
      }
    }
    res.add(prev);
    return res.toArray(new LongHit[res.size()]);
  }

  private LongHit[] getHits(String read, HashMap<String, AlnPos[]> index){
    int readLen = read.length();
    if(readLen < K){
      return new LongHit[0];
    }
    AlnPos[][] matches = new AlnPos[readLen - K + 1][];
    for(int i = 0; i <= readLen - K; i++){
      String s = read.substring(i, i + K);
      matches[i] = index.get(s);
      if(matches[i] == null){
        matches[i] = new AlnPos[0];
      }
    }
    return mergeHits(matches);
  }

  MappedRead[] getNBestRegions(byte[] read, HashMap<String, AlnPos[]> index, int from, int to, int minCount){
    String readStr = new String(read);
    LongHit[] hits = getHits(readStr, index);
    if(hits.length == 0){
      return null;
    }
    int readLen = read.length;
    int totalHit = hits.length;
    boolean[] isConcordant = concordanceArray(hits);

    Location location = getInitState(from, to, readLen, hits, isConcordant);
    ArrayList<MappedRead> bestPositions = new ArrayList<>();
    if(location.count > minCount){
      bestPositions.add(new MappedRead(location.start, location.end, location.count));
    }
    for(int i = 1; i < totalHit; i++){
      computeCount(from, to, readLen, hits, i, location, isConcordant);
      if(location.count > minCount){
        bestPositions.add(new MappedRead(location.start, location.end, location.count));
      }
    }
    if(bestPositions.size() == 0){
      return null;
    }

    return pruneWindows(bestPositions);
  }

  int computeKMerCount(byte[] read, HashMap<String, AlnPos[]> index, int from, int to){
    String readStr = new String(read);
    LongHit[] hits = getHits(readStr, index);
    if(hits.length == 0){
      return 0;
    }
    int readLen = read.length;
    int totalHit = hits.length;
    boolean[] isAdjacent = concordanceArray(hits);
    Location location = getInitState(from, to, readLen, hits, isAdjacent);
    int maxCount = 0;
    for(int i = 1; i < totalHit; i++){
      computeCount(from, to, readLen, hits, i, location, isAdjacent);
      if(location.count > maxCount){
        maxCount = location.count;
      }
    }
    return maxCount;
  }

}
