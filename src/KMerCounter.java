import org.jdom2.Document;
import org.jdom2.Element;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;

/**
 * Created by Gennady on 06.12.2015.
 */
class KMerCounter{

  int K;
  private String randomReadsCountsPath;
  private float pValue;
  int[] cutoffs;
  private float lowCoef;
  private float highCoef;
  private static KMerCounter instance;

  private LongHit[] mergeHits(int[][] hits){
    ArrayList<LongHit> longHits = new ArrayList();
    ArrayList<LongHit> currSet = new ArrayList<>();
    ArrayList<LongHit> nextSet = new ArrayList<>();
    int firstMatch = 0;
    do{
      if(firstMatch == hits.length){
        return new LongHit[0];
      }
      for(int j = 0; j < hits[firstMatch].length; j++){
        int pos_j = hits[firstMatch][j];
        LongHit longHit = new LongHit(pos_j, firstMatch);
        currSet.add(longHit);
        longHits.add(longHit);
      }
      firstMatch ++;
    }while(currSet.isEmpty());
    for(int i = firstMatch; i < hits.length; i++){
      int[] hits_i = hits[i];
      if(hits_i.length == 0){
        currSet = new ArrayList<>();
        nextSet = new ArrayList<>();
        continue;
      }
      if(currSet.size() == 0){
        for(int pos_j : hits_i){
          LongHit longHit = new LongHit(pos_j, i);
          currSet.add(longHit);
          longHits.add(longHit);
        }
        continue;
      }
      Iterator<LongHit> iteratorSet = currSet.iterator();
      LongHit lHit = iteratorSet.next();
      int j = 0;
      int rHit = hits_i[j];
      while(true){
        // no intersection
        if(lHit.genomePos + lHit.len + 1 < rHit){
          //longHits.add(lHit);
          if(iteratorSet.hasNext()){
            lHit = iteratorSet.next();
          }else{
            LongHit longHit = new LongHit(rHit, i);
            nextSet.add(longHit);
            longHits.add(longHit);
            j++;
            break;
          }
        }else if(lHit.genomePos + lHit.len + 1 > rHit){
          j++;
          LongHit longHit = new LongHit(rHit, i);
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
        LongHit longHit = new LongHit(rHit, i);
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

  private LongHit[] getHits(String read, HashMap<String, int[]> index){
    int readLen = read.length();
    if(readLen < K){
      return new LongHit[0];
    }
    int[][] matches = new int[readLen - K + 1][];
    for(int i = 0; i <= readLen - K; i++){
      String s = read.substring(i, i + K);
      matches[i] = index.get(s);
      if(matches[i] == null){
        matches[i] = new int[0];
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

  private Location getInitState(int[] contigEnds, int readLen,
                                LongHit[] longHits, boolean[] isAdjacent){
    Location res = new Location();
    int totalHit = longHits.length;
    LongHit h_i = longHits[0];
    int startMax = 0;
    int endMax = 0;
    for(int j = 0, start = 0; j < contigEnds.length; j++){
      int end = contigEnds[j];
      int locStart = Math.max(start, h_i.genomePos - h_i.readPos*3/2);
      int locEnd = Math.min(h_i.genomePos + h_i.len + K + (readLen - h_i.readPos - h_i.len - K)*3/2, end);
      if(locEnd - locStart > endMax - startMax){
        endMax = locEnd;
        startMax = locStart;
      }
      start = end;
    }
    res.start = startMax;
    res.end = endMax;
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

  private void computeCount(int[] contigEnds, int readLen, LongHit[] longHits,
                            int i, Location location, boolean[] isAdjacent){
    int totalHit = longHits.length;
    LongHit h_i = longHits[i];
    int start = 0;
    int end = 0;
    for(int j = 0, s = 0; j < contigEnds.length; j++){
      int e = contigEnds[j];
      int locStart = Math.max(s, h_i.genomePos - h_i.readPos*3/2);
      int locEnd = Math.min(h_i.genomePos + h_i.len + K + (readLen - h_i.readPos - h_i.len - K)*3/2, e);
      if(locEnd - locStart > end - start){
        end = locEnd;
        start = locStart;
      }
      s = e;
    }
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

  int computeKMerCount(byte[] read, HashMap<String, int[]> index, int from, int to){
    int readLen = read.length;
    String readStr = new String(read);
    LongHit[] matches = getHits(readStr, index);
    int totalHit = matches.length;
    if(totalHit == 0){
      return 0;
    }

    boolean[] isAdjacent = new boolean[totalHit];
    Hit prev = matches[0];
    for(int i = 1; i < totalHit; i++){
      Hit h_j = matches[i];
      float rate = (float)(h_j.genomePos - prev.genomePos)/(h_j.readPos - prev.readPos);
      if(lowCoef <= rate && rate <= highCoef){
        isAdjacent[i - 1] = true;
      }
      prev = h_j;
    }
    Location state = getInitState(from, to, readLen, matches, isAdjacent);
    int maxCount = state.count;
    for(int i = 1; i < totalHit; i++){
      computeCount(from, to, readLen, matches, i, state, isAdjacent);
      if(state.count > maxCount){
        maxCount = state.count;
      }
    }
    return maxCount;
  }


  private KMerCounter(Document document){
    try{
      Element element = document.getRootElement().getChild("Mapper");
      randomReadsCountsPath = element.getChildText("RandomModelPath");
      pValue = Float.parseFloat(element.getChildText("pValue"));
      cutoffs = readRandomModel();
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

  MappedRead mapReadToRegion(byte[] read, HashMap<String, int[]> index, int from, int to){
    int readLen = read.length;
    String readStr = new String(read);
    LongHit[] matches = getHits(readStr, index);
    int totalHit = matches.length;
    if(totalHit == 0){
      return null;
    }
    boolean[] isAdjacent = new boolean[totalHit];
    LongHit prev = matches[0];
    for(int i = 1; i < totalHit; i++){
      LongHit h_j = matches[i];
      float rate = (float)(h_j.genomePos - prev.genomePos)/(h_j.readPos - prev.readPos);
      if(lowCoef <= rate && rate <= highCoef){
        isAdjacent[i - 1] = true;
      }
      prev = h_j;
    }
    Location state = getInitState(from, to, readLen, matches, isAdjacent);
    MappedRead res = new MappedRead(state.start, state.end, state.count);
    for(int i = 1; i < totalHit; i++){
      computeCount(from, to, readLen, matches, i, state, isAdjacent);
      if(state.count > res.count){
        res.count = state.count;
        res.start = state.start;
        res.end = state.end;
      }
    }
    res.seq = read;
    return res;
  }

  MappedRead mapReadToRegion(byte[] read, HashMap<String, int[]> index, int[] contigEnds){
    int readLen = read.length;
    String readStr = new String(read);
    LongHit[] matches = getHits(readStr, index);
    int totalHit = matches.length;
    if(totalHit == 0){
      return null;
    }
    boolean[] isAdjacent = new boolean[totalHit];
    Hit prev = matches[0];
    for(int i = 1; i < totalHit; i++){
      Hit h_j = matches[i];
      float rate = (float)(h_j.genomePos - prev.genomePos)/(h_j.readPos - prev.readPos);
      if(lowCoef <= rate && rate <= highCoef){
        isAdjacent[i - 1] = true;
      }
      prev = h_j;
    }
    Location state = getInitState(contigEnds, readLen, matches, isAdjacent);
    MappedRead res = new MappedRead(state.start, state.end, state.count);
    for(int i = 1; i < totalHit; i++){
      computeCount(contigEnds, readLen, matches, i, state, isAdjacent);
      if(state.count > res.count){
        res.count = state.count;
        res.start = state.start;
        res.end = state.end;
      }
    }
    res.seq = read;
    return res;
  }

  MappedRead[] getNBestRegions(byte[] read, HashMap<String, int[]> index, int from, int to){
    int readLen = read.length;
    String readStr = new String(read);
    LongHit[] matches = getHits(readStr, index);
    int totalHit = matches.length;
    if(totalHit == 0){
      return null;
    }
    boolean[] isAdjacent = new boolean[totalHit];
    Hit prev = matches[0];
    for(int i = 1; i < totalHit; i++){
      Hit h_j = matches[i];
      float rate = (float)(h_j.genomePos - prev.genomePos)/(h_j.readPos - prev.readPos);
      if(lowCoef <= rate && rate <= highCoef){
        isAdjacent[i - 1] = true;
      }
      prev = h_j;
    }
    int cutoff = cutoffs[readLen];
    Location location = getInitState(from, to, readLen, matches, isAdjacent);
    ArrayList<MappedRead> bestPositions = new ArrayList<>();
    if(location.count > cutoff){
      bestPositions.add(new MappedRead(location.start, location.end, location.count));
    }
    for(int i = 1; i < totalHit; i++){
      computeCount(from, to, readLen, matches, i, location, isAdjacent);
      if(location.count > cutoff){
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

  MappedRead[] getNBestRegions(byte[] read, HashMap<String, int[]> index, int[] contigEnds){
    int readLen = read.length;
    String readStr = new String(read);
    LongHit[] matches = getHits(readStr, index);
    int totalHit = matches.length;
    if(totalHit == 0){
      return null;
    }
    boolean[] isAdjacent = new boolean[totalHit];
    Hit prev = matches[0];
    for(int i = 1; i < totalHit; i++){
      Hit h_j = matches[i];
      float rate = (float)(h_j.genomePos - prev.genomePos)/(h_j.readPos - prev.readPos);
      if(lowCoef <= rate && rate <= highCoef){
        isAdjacent[i - 1] = true;
      }
      prev = h_j;
    }
    int cutoff = cutoffs[readLen];
    Location location = getInitState(contigEnds, readLen, matches, isAdjacent);
    ArrayList<MappedRead> bestPositions = new ArrayList<>();
    if(location.count > cutoff){
      bestPositions.add(new MappedRead(location.start, location.end, location.count));
    }
    for(int i = 1; i < totalHit; i++){
      computeCount(contigEnds, readLen, matches, i, location, isAdjacent);
      if(location.count > cutoff){
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

}
