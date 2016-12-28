import org.jdom2.Document;
import org.jdom2.Element;

import java.io.*;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * Date: 03.01.15
 */
class KMerCounterMSA extends Constants{

  int K;
  private String randomReadsCountsPath;
  private float pValue;
  int[] cutoffs;
  private float lowCoef;
  private float highCoef;
  private static KMerCounterMSA instance;

  private int[] readRandomModel()
      throws IOException{
    DataInputStream inputStream = new DataInputStream(new FileInputStream(randomReadsCountsPath));
    K = inputStream.readInt();
    highCoef = inputStream.readFloat();
    lowCoef = 1/highCoef;
    int minReadLen = inputStream.readInt();
    int maxReadLen = inputStream.readInt();
    int step = inputStream.readInt();
    int num = inputStream.readInt();
    int size = (maxReadLen - minReadLen)/step + 1;
    int[] cutoffsInterpolation = new int[size];
    for(int i = 0; i < size; i++){
      int[] counts = new int[num];
      for(int j = 0; j < num; j++){
        counts[j] = inputStream.readInt();
      }
      int index = Math.round(num*(1.0f - pValue));
      cutoffsInterpolation[i] = counts[index];
    }
    inputStream.close();
    int[] cutoffs = new int[maxReadLen + 1];
    for(int i = 0; i < minReadLen; i++){
      cutoffs[i] = cutoffsInterpolation[0];
    }
    for(int i = minReadLen; i < maxReadLen; i++){
      int left = (i - minReadLen)/step;
      int distToLeft = i - left*step - minReadLen;
      int distToRight = step - distToLeft;
      cutoffs[i] = (distToRight*cutoffsInterpolation[left] +
          distToLeft*cutoffsInterpolation[left + 1])/step;
    }
    cutoffs[maxReadLen] = cutoffsInterpolation[size - 1];
    return cutoffs;
  }

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

  private class CoordComparator implements Comparator{

    @Override
    public int compare(Object o1, Object o2){
      MappedRead mRead1 = (MappedRead) o1;
      MappedRead mRead2 = (MappedRead) o2;
      if(mRead1.start == mRead2.start){
        return mRead2.end - mRead1.end;
      }
      return mRead1.start - mRead2.start;
    }
  }

  private LongHit[] mergeHits(AlnPos[][] hits){
    ArrayList<LongHit> longHits = new ArrayList();
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

  private Location getInitState(int from, int to, int readLen,
                                LongHit[] longHits, boolean[] isAdjacent){
    Location res = new Location();
    int totalHit = longHits.length;
    LongHit h_i = longHits[0];
    res.start = Math.max(from, h_i.genomePos - h_i.readPos*3/2);
    res.end = Math.min(h_i.genomePos + h_i.len + K + (readLen - h_i.readPos - h_i.len - K)*3/2, to);
    res.startIndex = 0;
    int count = h_i.len;
    int j;
    for(j = 1; j < totalHit ; j++){
      LongHit h_j = longHits[j];
      if(h_j.genomePos + K + h_j.len >= res.end){
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

  private void computeCount(int from, int to, int readLen, LongHit[] longHits,
                            int i, Location location, boolean[] isAdjacent){
    int totalHit = longHits.length;
    LongHit h_i = longHits[i];
    int start = Math.max(from, h_i.genomePos - h_i.readPos*3/2);
    int end = Math.min(h_i.genomePos + h_i.len + K + (readLen - h_i.readPos - h_i.len - K)*3/2, to);
    int j;
    int count = location.count;
    if(start != location.start){
      int sIndex;
      if(start < location.start){
        for(j = location.startIndex - 1; j >= 0; j--){
          LongHit h_j = longHits[j];
          if(h_j.genomePos < start){
            break;
          }
          if(isAdjacent[j]){
            count ++;
          }
          count += h_j.len;
        }
        sIndex = j + 1;
      }else{
        LongHit h_j = longHits[location.startIndex];
        sIndex = location.startIndex;
        if(h_j.genomePos < start){
          count -= h_j.len;
          for(j = location.startIndex + 1; j < totalHit; j++){
            h_j = longHits[j];
            if(isAdjacent[j - 1]){
              count --;
            }
            if(h_j.genomePos >= start){
              break;
            }
            count -= h_j.len;
          }
          if(j < totalHit){
            sIndex = j;
          }
        }
      }
      location.start = start;
      location.startIndex = sIndex;
    }
    if(end != location.end){
      int eIndex;
      if(end > location.end){
        for(j = location.endIndex + 1; j < totalHit ; j++){
          LongHit h_j = longHits[j];
          if(h_j.genomePos + K + h_j.len > end){
            break;
          }
          if(isAdjacent[j - 1]){
            count ++;
          }
          count += h_j.len;
        }
        eIndex = j - 1;
      }else{
        LongHit h_j = longHits[location.endIndex];
        eIndex = location.endIndex;
        if(h_j.genomePos + K + h_j.len > end){
          count -= h_j.len;
          for(j = location.endIndex - 1; j >= 0; j--){
            h_j = longHits[j];
            if(isAdjacent[j]){
              count --;
            }
            if(h_j.genomePos + K + h_j.len <= end){
              break;
            }
            count -= h_j.len;
          }
          if(j >= 0){
            eIndex = j;
          }
        }
      }
      location.end = end;
      location.endIndex = eIndex;
    }
    location.count = count;
  }

  MappedRead mapReadToRegion(byte[] read, HashMap<String, AlnPos[]> index, int from, int to){
    String readStr = new String(read);
    LongHit[] hits = getHits(readStr, index);
    if(hits.length == 0){
      return null;
    }
    int readLen = read.length;
    int totalHit = hits.length;
    boolean[] isAdjacent = new boolean[totalHit];
    LongHit prev = hits[0];
    for(int i = 1; i < totalHit; i++){
      LongHit h_j = hits[i];
      float rate = (float)(h_j.genomePos - prev.genomePos)/(h_j.readPos - prev.readPos);
      if(lowCoef <= rate && rate <= highCoef){
        isAdjacent[i - 1] = true;
      }
      prev = h_j;
    }
    Location location = getInitState(from, to, readLen, hits, isAdjacent);
    Location locationMax = new Location(location);
    for(int i = 1; i < totalHit; i++){
      computeCount(from, to, readLen, hits, i, location, isAdjacent);
      if(location.count > locationMax.count){
        locationMax = new Location(location);
      }
    }
    MappedRead res = new MappedRead(locationMax.start, locationMax.end, locationMax.count);
//    if(storeHits){
//      res.hits = new LongHit[locationMax.endIndex - locationMax.startIndex];
//      System.arraycopy(hits, locationMax.startIndex, res.hits, 0, res.hits.length);
//    }
    res.seq = read;
    return res;
  }

  MappedRead[] getNBestRegions(byte[] read, HashMap<String, AlnPos[]> index, int from, int to, int minCount){
    String readStr = new String(read);
    LongHit[] hits = getHits(readStr, index);
    if(hits.length == 0){
      return null;
    }
    int readLen = read.length;
    int totalHit = hits.length;
    boolean[] isAdjacent = new boolean[totalHit];
    LongHit prev = hits[0];
    for(int i = 1; i < totalHit; i++){
      LongHit h_j = hits[i];
      float rate = (float)(h_j.genomePos - prev.genomePos)/(h_j.readPos - prev.readPos);
      if(lowCoef <= rate && rate <= highCoef){
        isAdjacent[i - 1] = true;
      }
      prev = h_j;
    }
    Location location = getInitState(from, to, readLen, hits, isAdjacent);
    ArrayList<MappedRead> bestPositions = new ArrayList<>();
    if(location.count > minCount){
      bestPositions.add(new MappedRead(location.start, location.end, location.count));
    }
    for(int i = 1; i < totalHit; i++){
      computeCount(from, to, readLen, hits, i, location, isAdjacent);
      if(location.count > minCount){
        bestPositions.add(new MappedRead(location.start, location.end, location.count));
      }
    }
    if(bestPositions.size() == 0){
      return null;
    }
    Collections.sort(bestPositions, new CoordComparator());
    ArrayList<MappedRead> res = new ArrayList<>();
    MappedRead pivot = bestPositions.get(0);
    for(int i = 1; i < bestPositions.size(); i++){
      MappedRead curr = bestPositions.get(i);
      if(pivot.end > curr.start){
        // intersection
        if(pivot.count < curr.count){
          pivot = curr;
        }
      }else{
        res.add(pivot);
        pivot = curr;
      }
    }
    res.add(pivot);
    return res.toArray(new MappedRead[res.size()]);
  }

  int computeKMerCount(byte[] read, HashMap<String, AlnPos[]> index, int from, int to){
    String readStr = new String(read);
    LongHit[] hits = getHits(readStr, index);
    if(hits.length == 0){
      return 0;
    }
    int readLen = read.length;
    int totalHit = hits.length;
    boolean[] isAdjacent = new boolean[totalHit];
    LongHit prev = hits[0];
    for(int i = 1; i < totalHit; i++){
      LongHit h_j = hits[i];
      float rate = (float)(h_j.genomePos - prev.genomePos)/(h_j.readPos - prev.readPos);
      if(lowCoef <= rate && rate <= highCoef){
        isAdjacent[i - 1] = true;
      }
      prev = h_j;
    }
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
