import java.util.HashMap;

/**
 * Created by Геннадий on 06.08.2014.
 */
class Interval extends Constants implements Comparable{

  int start;
  int end;
  short junction;
  byte[] concat;
  byte[] left;
  byte[] right;
  String seqFinal;
  int[] leftCounts;
  int[] rightCounts;
  boolean merge;
  boolean updatedLeft;
  boolean updatedRight;
  HashMap<String, int[]> index;
  HashMap<String, int[]> indexLeft;
  HashMap<String, int[]> indexRight;
  boolean updateCount;

  int finalizingReadsFound;

  void init(int K, int threadNum) throws InterruptedException{
    leftCounts = new int[4];
    rightCounts = new int[4];
    index = new HashMap<>();
    indexLeft = new HashMap<>();
    indexRight = new HashMap<>();
    buildIndexConcat(K, threadNum);
    buildIndexLeft(K, threadNum);
    buildIndexRight(K, threadNum);
  }

  @Override
  public int compareTo(Object o){
    Interval interval = (Interval) o;
    return start - interval.start;
  }

  boolean growLeftEdge(int coverageThreshold){
    if(seqFinal != null){
      return false;
    }
    updatedLeft = false;
    int maxCount = leftCounts[0];
    int maxIndex = 0;
    leftCounts[0] = 0;
    for(int i = 1; i < 4; i++){
      if(leftCounts[i] > maxCount){
        maxCount = leftCounts[i];
        maxIndex = i;
      }
      leftCounts[i] = 0;
    }
    if(maxCount > coverageThreshold){
      byte[] concatN = new byte[concat.length + 1];
      byte[] leftN = new byte[left.length + 1];
      System.arraycopy(left, 0, concatN, 0, left.length);
      concatN[junction + 1] = iToN[maxIndex];
      System.arraycopy(right, 0, concatN, junction + 2, right.length);
      System.arraycopy(left, 0, leftN, 0, left.length);
      leftN[left.length] = iToN[maxIndex];
      concat = concatN;
      left = leftN;
      junction++;
      updatedLeft = true;
      finalizingReadsFound = 0;
      return true;
    }
    return false;
  }

  boolean growRightEdge(int coverageThreshold){
    if(seqFinal != null){
      return false;
    }
    updatedRight = false;
    int maxCount = rightCounts[0];
    int maxIndex = 0;
    rightCounts[0] = 0;
    for(int i = 1; i < 4; i++){
      if(rightCounts[i] > maxCount){
        maxCount = rightCounts[i];
        maxIndex = i;
      }
      rightCounts[i] = 0;
    }
    if(maxCount > coverageThreshold){
      byte[] concatN = new byte[concat.length + 1];
      byte[] rightN = new byte[right.length + 1];
      System.arraycopy(left, 0, concatN, 0, left.length);
      concatN[junction + 1] = iToN[maxIndex];
      System.arraycopy(right, 0, concatN, junction + 2,
          right.length);
      System.arraycopy(right, 0, rightN, 1, right.length);
      rightN[0] = iToN[maxIndex];
      concat = concatN;
      right = rightN;
      updatedRight = true;
      finalizingReadsFound = 0;
      return true;
    }
    return false;
  }

  void buildIndexLeft(int K, int threadNum) throws InterruptedException{
    indexLeft = StringIndexer.buildIndexPara(new String(left), K, threadNum);
  }

  void buildIndexRight(int K, int threadNum) throws InterruptedException{
    indexRight = StringIndexer.buildIndexPara(new String(right), K, threadNum);
  }

  void buildIndexConcat(int K, int threadNum) throws InterruptedException{
    index = StringIndexer.buildIndexPara(new String(concat), K, threadNum);
  }

}
