import gnu.trove.iterator.TIntFloatIterator;
import gnu.trove.map.hash.TIntFloatHashMap;
import org.jdom2.Document;
import org.jdom2.Element;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Created with IntelliJ IDEA.
 * Date: 08.08.15
 */
class MinRefSetSelector extends Constants{

  private float delta;
  private int minContigLen;
  private int match;
  private int mismatch;
  private int gop;//gap open penalty
  private int gep; //gap extension penalty
  private ReferenceAlignment refAlignment;
  private String outPath;
  private Logger logger;
  private int threadNum;
  private int maxNongreedyComponentNum;

  MinRefSetSelector(Document document){
    Element element = document.getRootElement().getChild("ReferenceSelector");
    outPath = document.getRootElement().getChildText("OutPath");
    minContigLen = Integer.parseInt(element.getChildText("MinContigLength"));
    delta = Float.parseFloat(element.getChildText("Delta"));
    Element alignerElement = document.getRootElement().getChild("Mapper").getChild("Aligner");
    match = Integer.parseInt(alignerElement.getChildText("Match"));
    mismatch = Integer.parseInt(alignerElement.getChildText("Mismatch"));
    gop = Integer.parseInt(alignerElement.getChildText("GapOpenPenalty"));
    gep = Integer.parseInt(alignerElement.getChildText("GapExtensionPenalty"));
    maxNongreedyComponentNum = Integer.parseInt(element.getChildText("MaxNongreedyComponentNumber"));
    refAlignment = ReferenceAlignment.getInstance(document);
    logger = Logger.getInstance(document);
    threadNum = Integer.parseInt(document.getRootElement().getChildText("ThreadNumber"));
    if(threadNum == -1){
      threadNum = Runtime.getRuntime().availableProcessors();
    }
  }

  private class Score implements Comparable{
    public int score;
    public float identity;
    public Reference ref;

    @Override
    public int compareTo(Object o){
      Score sc = (Score) o;
      if(sc.identity == identity){
        return 0;
      }
      return (sc.identity > identity)?1:-1;
    }

    public String toString(){
      StringBuilder builder = new StringBuilder();
      builder.append(ref.name).append(' ').append(identity).append('\n');
      return builder.toString();
    }
  }

  private class ReferenceSelector implements Runnable{

    public Score[] scores;
    public byte[] contig;
    public int contigAlnStart;
    public int contigAlnEnd;
    public Path path;

    @Override
    public void run(){
      scores = new Score[refAlignment.refAlns.size()];
      for(int k = 0; k < scores.length; k++){
        Reference seq = refAlignment.refAlns.get(k);
        int score = 0;
        int identity = 0;
        int len = 0;
        boolean gapOpened = false;
        for(int i = contigAlnStart, j = 0; i < contigAlnEnd; i++, j++){
          if(contig[j] == GAP){
            if(seq.aln.charAt(i) != GAP){
              if(gapOpened){
                score -= gep;
              }else{
                score -= gop;
                gapOpened = true;
              }
              len ++;
            }
          }else{
            if(seq.aln.charAt(i) == GAP){
              if(gapOpened){
                score -= gep;
              }else{
                score -= gop;
                gapOpened = true;
              }
              len ++;
            }else{
              gapOpened = false;
              if(contig[j] == seq.aln.charAt(i)){
                score += match;
                identity ++;
              }else{
                score += mismatch;
              }
              len ++;
            }
          }
        }
        Score sc = new Score();
        sc.score = score;
        sc.ref = seq;
        sc.identity = (float)identity/len;
        scores[k] = sc;
      }
      Arrays.sort(scores);
    }
  }

  private HashMap<String, TIntFloatHashMap> selectRefRecursive(
      HashMap<String, TIntFloatHashMap> refToContigIds, int maxCount, ArrayList<Path> longestPaths){
    HashMap<String, TIntFloatHashMap> selectedReferences = null;
    int minCount = Integer.MAX_VALUE;
    float maxSum = 0;
    for(Map.Entry<String, TIntFloatHashMap> entry : refToContigIds.entrySet()){
      String ref = entry.getKey();
      TIntFloatHashMap map = entry.getValue();
      HashMap<String, TIntFloatHashMap> selReferences = new HashMap<>();
      selReferences.put(ref, map);

      HashMap<String, TIntFloatHashMap> refToContigIdsRest = new HashMap<>();
      for(Map.Entry<String, TIntFloatHashMap> ent: refToContigIds.entrySet()){
        String r = ent.getKey();
        if(r.equals(ref)){
          continue;
        }
        TIntFloatHashMap m = ent.getValue();
        TIntFloatHashMap mRest = new TIntFloatHashMap();
        for(TIntFloatIterator iter = m.iterator(); iter.hasNext(); ){
          iter.advance();
          if(!map.contains(iter.key())){
            mRest.put(iter.key(), iter.value());
          }
        }
        if(!mRest.isEmpty()){
          refToContigIdsRest.put(r, mRest);
        }
      }
      if(maxCount > 1){
        if(!refToContigIdsRest.isEmpty()){
          HashMap<String, TIntFloatHashMap> selRef = selectRefRecursive(refToContigIdsRest, maxCount - 1, longestPaths);
          if(selRef != null){
            selReferences.putAll(selRef);
          }else{
            continue;
          }
        }
      }else{
        if(!refToContigIdsRest.isEmpty()){
          continue;
        }
      }

      int count = selReferences.size();
      float sum = 0;
      int totalLen = 0;
      for(TIntFloatHashMap map1: selReferences.values()){
        for(TIntFloatIterator iter = map1.iterator(); iter.hasNext();){
          iter.advance();
          int len = longestPaths.get(iter.key()).contig.length;
          sum += iter.value()*len;
          totalLen += len;
        }
      }
      sum /= totalLen;
      if(count < minCount){
        minCount = count;
        maxSum = sum;
        selectedReferences = selReferences;
      }else if(count == minCount){
        if(sum > maxSum){
          maxSum = sum;
          selectedReferences = selReferences;
        }
      }
    }
    return selectedReferences;
  }

  private HashMap<String, int[]> selectRefSubset(ArrayList<ReferenceSelector> tasks, ArrayList<Path> longestPaths){
    HashMap<String, TIntFloatHashMap> refToContigIds = new HashMap<>();
    for(ReferenceSelector task : tasks){
      Score[] scores = task.scores;
      float maxIdentity = scores[0].identity;
      for(Score score : scores){
        if(maxIdentity - score.identity < delta){
          TIntFloatHashMap map = refToContigIds.get(score.ref.name);
          if(map == null){
            map = new TIntFloatHashMap();
            refToContigIds.put(score.ref.name, map);
          }
          map.put(task.path.id, score.identity);
        }
      }
    }
    HashMap<String, int[]> selectedReferences = selectReferencesGreedy(refToContigIds, longestPaths);
    if(selectedReferences.size() > maxNongreedyComponentNum){
      return selectedReferences;
    }
    if(selectedReferences.size() > 1){
      HashMap<String, TIntFloatHashMap> selRef =
          selectRefRecursive(refToContigIds, selectedReferences.size(), longestPaths);
      selectedReferences.clear();
      for(Map.Entry<String, TIntFloatHashMap> entry : selRef.entrySet()){
        selectedReferences.put(entry.getKey(), entry.getValue().keys());
      }
    }
    return selectedReferences;
  }

  private HashMap<String, int[]> selectReferencesGreedy(
      HashMap<String, TIntFloatHashMap> refToContigIds, ArrayList<Path> longestPaths){
    HashMap<String, int[]> selectedReferences = new HashMap<>();
    HashMap<String, TIntFloatHashMap> rToContigIds = new HashMap<>();
    for(Map.Entry<String, TIntFloatHashMap> entry: refToContigIds.entrySet()){
      TIntFloatHashMap map = new TIntFloatHashMap(entry.getValue());
      rToContigIds.put(entry.getKey(), map);
    }
    while(!rToContigIds.isEmpty()){
      int maxCount = 0;
      float maxSum = 0;
      String maxRef = null;
      for(Map.Entry<String, TIntFloatHashMap> entry: rToContigIds.entrySet()){
        TIntFloatHashMap map = entry.getValue();
        if(map.size() > maxCount){
          maxCount = map.size();
          maxRef = entry.getKey();
          maxSum = 0;
          int totalLen = 0;
          for(TIntFloatIterator iter = map.iterator();iter.hasNext();){
            iter.advance();
            int len = longestPaths.get(iter.key()).contig.length;
            maxSum += iter.value()*len;
            totalLen += len;
          }
          maxSum /= totalLen;
        }else if(map.size() == maxCount){
          float currSum = 0;
          int totalLen = 0;
          for(TIntFloatIterator iter = map.iterator();iter.hasNext();){
            iter.advance();
            int len = longestPaths.get(iter.key()).contig.length;
            currSum += iter.value()*len;
            totalLen += len;
          }
          currSum /= totalLen;
          if(currSum > maxSum){
            maxSum = currSum;
            maxRef = entry.getKey();
          }
        }
      }
      TIntFloatHashMap maxMap = rToContigIds.get(maxRef);
      selectedReferences.put(maxRef, maxMap.keys());
      rToContigIds.remove(maxRef);
      ArrayList<String> del = new ArrayList<>();
      for(Map.Entry<String, TIntFloatHashMap> entry: rToContigIds.entrySet()){
        TIntFloatHashMap map = entry.getValue();
        for(TIntFloatIterator iter = maxMap.iterator();iter.hasNext();){
          iter.advance();
          map.remove(iter.key());
        }
        if(map.isEmpty()){
          del.add(entry.getKey());
        }
      }
      for(String ref: del){
        rToContigIds.remove(ref);
      }
    }
    return selectedReferences;
  }

  HashMap<String, int[]> chooseReferences(ArrayList<Path> longestPaths) throws InterruptedException, IOException{
    ArrayList<ReferenceSelector> tasks = new ArrayList<>();
    ExecutorService executor = Executors.newFixedThreadPool(threadNum);
    for(Path path: longestPaths){
      if(path.contig.length >= minContigLen){
        ReferenceSelector task = new ReferenceSelector();
        task.contig = path.contigAln;
        task.contigAlnStart = path.contigAlnStart;
        task.contigAlnEnd = path.contigAlnEnd;
        task.path = path;
        executor.execute(task);
        tasks.add(task);
      }
    }
    executor.shutdown();
    executor.awaitTermination(10, TimeUnit.HOURS);
    BufferedWriter writer = new BufferedWriter(new FileWriter(outPath + "references.txt"));
    for(ReferenceSelector task : tasks){
      writer.write("Path " + Integer.toString(task.path.id) + "\n");
      for(Score score : task.scores){
        writer.write(score.toString());
      }
    }
    writer.close();
    HashMap<String, int[]> selectedRef = selectRefSubset(tasks, longestPaths);
    logger.println("Selected references (with path IDs):");
    for(Map.Entry<String, int[]> entry: selectedRef.entrySet()){
      logger.printf(entry.getKey());
      for(int pathID: entry.getValue()){
        logger.printf(" %d", pathID);
      }
      logger.printf("\n");
    }
    return selectedRef;
  }


}
