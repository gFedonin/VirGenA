import java.util.*;

class KMerCounterBase extends Constants{

  int K;
//  String randomReadsCountsPath;
  float pValue;
  int[] cutoffs;
  float lowCoef;
  float highCoef;


  Location getInitState(int from, int to, int readLen,
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

  private void updateCount(int start, int end, LongHit[] longHits,
                   Location location, boolean[] isAdjacent){
    int totalHit = longHits.length;
    int j;
    int count = location.count;
    int sIndex;
    if(start != location.start){
      sIndex = -1;
      boolean prevInWindow;
      LongHit startHit = longHits[location.startIndex];
      if(start < location.start){
        prevInWindow = startHit.genomePos + K + startHit.len <= end;
        for(j = location.startIndex - 1; j >= 0; j--){
          LongHit h_j = longHits[j];
          if(h_j.genomePos < start){
            sIndex = j + 1;
            break;
          }
          if(h_j.genomePos + K + h_j.len > end){
            continue;
          }
          if(prevInWindow && isAdjacent[j]){
            count ++;
          }
          count += h_j.len;
          prevInWindow = true;
        }
        if(sIndex == -1){
          sIndex = 0;
        }
      }else{
        LongHit h_j = longHits[location.startIndex];
        prevInWindow = h_j.genomePos >= start;
        if(prevInWindow){
          sIndex = location.startIndex;
        }else{
          count -= h_j.len;
          for(j = location.startIndex + 1; j <= location.endIndex; j++){
            h_j = longHits[j];
            if(isAdjacent[j - 1]){
              count --;
            }
            if(h_j.genomePos >= start){
              break;
            }
            count -= h_j.len;
          }
          if(j <= location.endIndex){
            sIndex = j;
          }else{
            for(sIndex = j; sIndex < totalHit; sIndex++){
              h_j = longHits[sIndex];
              if(h_j.genomePos >= start){
                break;
              }
            }
          }
        }
      }
      location.start = start;
    }else{
      sIndex = location.startIndex;
    }
    if(end != location.end){
      int eIndex = -1;
      boolean prevInWindow;
      if(end > location.end){
        LongHit endHit = longHits[location.endIndex];
        prevInWindow = endHit.genomePos >= start;
        for(j = location.endIndex + 1; j < totalHit ; j++){
          LongHit h_j = longHits[j];
          if(h_j.genomePos + K + h_j.len > end){
            eIndex = j - 1;
            break;
          }
          if(h_j.genomePos < start){
            continue;
          }
          if(prevInWindow && isAdjacent[j - 1]){
            count ++;
          }
          count += h_j.len;
          prevInWindow = true;
        }
        if(eIndex == -1){
          eIndex = totalHit - 1;
        }
      }else{
        LongHit h_j = longHits[location.endIndex];
        if(h_j.genomePos + K + h_j.len > end){
          count -= h_j.len;
          for(j = location.endIndex - 1; j >= location.startIndex; j--){
            h_j = longHits[j];
            if(isAdjacent[j]){
              count --;
            }
            if(h_j.genomePos + K + h_j.len <= end){
              break;
            }
            count -= h_j.len;
          }
          if(j >= location.startIndex){
            eIndex = j;
          }else{
            for(eIndex = j; eIndex >= 0; eIndex--){
              h_j = longHits[eIndex];
              if(h_j.genomePos + K + h_j.len <= end){
                break;
              }
            }
          }
        }else{
          eIndex = location.endIndex;
        }
      }
      location.end = end;
      location.endIndex = eIndex;
    }
    location.startIndex = sIndex;
    location.count = count;
  }

  void computeCount(int from, int to, int readLen, LongHit[] longHits,
                    int i, Location location, boolean[] isAdjacent){
    LongHit h_i = longHits[i];
    int start = Math.max(from, h_i.genomePos - h_i.readPos*3/2);
    int end = Math.min(h_i.genomePos + h_i.len + K + (readLen - h_i.readPos - h_i.len - K)*3/2, to);
    updateCount(start, end, longHits, location, isAdjacent);
  }


  void computeCount(int[] contigEnds, int readLen, LongHit[] longHits,
                            int i, Location location, boolean[] isAdjacent){
    LongHit h_i = longHits[i];
    int start = 0;
    int end = 0;
    for(int j = 0, s = 0; j < contigEnds.length; j++){
      int e = contigEnds[j];
      if(s <= h_i.genomePos && h_i.genomePos < e){
        start = Math.max(s, h_i.genomePos - h_i.readPos*3/2);
        end = Math.min(h_i.genomePos + h_i.len + K + (readLen - h_i.readPos - h_i.len - K)*3/2, e);
        if(h_i.genomePos + h_i.len + K > e){
          return;
        }
        break;
      }
//      if(locEnd - locStart > end - start){
//        end = locEnd;
//        start = locStart;
//      }
      s = e;
    }
    updateCount(start, end, longHits, location, isAdjacent);
  }


//  void readRandomModel()
//      throws IOException{
//    DataInputStream inputStream = new DataInputStream(new FileInputStream(randomReadsCountsPath));
//    K = inputStream.readInt();
//    highCoef = inputStream.readFloat();
//    lowCoef = 1/highCoef;
//    int minReadLen = inputStream.readInt();
//    int maxReadLen = inputStream.readInt();
//    int step = inputStream.readInt();
//    int num = inputStream.readInt();
//    int size = (maxReadLen - minReadLen)/step + 1;
//    int[] cutoffsInterpolation = new int[size];
//    for(int i = 0; i < size; i++){
//      int[] counts = new int[num];
//      for(int j = 0; j < num; j++){
//        counts[j] = inputStream.readInt();
//      }
//      int index = Math.round(num*(1.0f - pValue));
//      cutoffsInterpolation[i] = counts[index];
//    }
//    inputStream.close();
//    cutoffs = new int[maxReadLen + 1];
//    for(int i = 0; i < minReadLen; i++){
//      cutoffs[i] = cutoffsInterpolation[0];
//    }
//    for(int i = minReadLen; i < maxReadLen; i++){
//      int left = (i - minReadLen)/step;
//      int distToLeft = i - left*step - minReadLen;
//      int distToRight = step - distToLeft;
//      cutoffs[i] = (distToRight*cutoffsInterpolation[left] +
//          distToLeft*cutoffsInterpolation[left + 1])/step;
//    }
//    cutoffs[maxReadLen] = cutoffsInterpolation[size - 1];
//  }

  void computeCuttofs(int minReadLen, int maxReadLen, int step, int[][] counts){
    int num = counts.length;
    int size = (maxReadLen - minReadLen)/step + 1;
    int[] cutoffsInterpolation = new int[size];
    for(int i = 0; i < size; i++){
      int index = Math.round(num*(1.0f - pValue));
      cutoffsInterpolation[i] = counts[i][index];
    }
    cutoffs = new int[maxReadLen + 1];
    for(int i = 0; i <= minReadLen; i++){
      cutoffs[i] = cutoffsInterpolation[0];
    }
    for(int i = maxReadLen; i > minReadLen ; i--){
      int right = size - (maxReadLen - i)/step - 1;
      if(right == 0){
        cutoffs[i] = cutoffsInterpolation[0];
      }else{
        int distToRight = (maxReadLen - i)%step;
        int distToLeft = step - distToRight;
        cutoffs[i] = (distToRight*cutoffsInterpolation[right - 1] + distToLeft*cutoffsInterpolation[right])/step;
      }
    }
  }

  class CoordComparator implements Comparator{

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

  boolean[] concordanceArray(LongHit[] longHits){
    int totalHit = longHits.length;
    boolean[] isAdjacent = new boolean[totalHit];
    Hit prev = longHits[0];
    for(int i = 1; i < totalHit; i++){
      Hit h_j = longHits[i];
      float rate = (float)(h_j.genomePos - prev.genomePos)/(h_j.readPos - prev.readPos);
      if(lowCoef <= rate && rate <= highCoef){
        isAdjacent[i - 1] = true;
      }
      prev = h_j;
    }
    return isAdjacent;
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

  MappedRead[] pruneWindows(ArrayList<MappedRead> bestPositions){
    bestPositions.sort(new CoordComparator());
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

  private LongHit[] mergeHits(int[][] hits){
    ArrayList<LongHit> longHits = new ArrayList<>();
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

  LongHit[] getHits(String read, HashMap<String, int[]> index){
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

  int computeKMerCount(byte[] read, HashMap<String, int[]> index, int from, int to){
    int readLen = read.length;
    String readStr = new String(read);
    LongHit[] matches = getHits(readStr, index);
    int totalHit = matches.length;
    if(totalHit == 0){
      return 0;
    }

    boolean[] isAdjacent = concordanceArray(matches);
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

}
