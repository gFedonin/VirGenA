import gnu.trove.list.array.TIntArrayList;

import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Created by Геннадий on 08.12.2014.
 */
class StringIndexer extends Constants{

  private static class StringIndexing implements Runnable{

    private int from;
    private int to;
    public HashMap<String, TIntArrayList> index;
    private String sequence;
    private int K;

    public StringIndexing(int from, int to, String sequence, int K){
      this.from = from;
      this.to = to;
      this.sequence = sequence;
      this.K = K;
    }

    @Override
    public void run(){
      index = new HashMap<>();
      if(to > sequence.length() - K){
        to = sequence.length() - K;
      }
      for(int i = from; i <= to; i++){
        String s = sequence.substring(i, i + K);
        TIntArrayList l = index.get(s);
        if(l == null){
          l = new TIntArrayList();
          index.put(s, l);
        }
        l.add(i);
      }
    }
  }

  static HashMap<String, int[]> buildIndexPara(String s, int K, int threadNum) throws InterruptedException{
    StringIndexing[] tasks = new StringIndexing[threadNum];
    ExecutorService executor = Executors.newFixedThreadPool(threadNum);
    int size = s.length()/threadNum;
    for(int i = 0; i < threadNum - 1; i++){
      StringIndexing task = new StringIndexing(i*size, (i + 1)*size, s, K);
      executor.execute(task);
      tasks[i] = task;
    }
    StringIndexing task = new StringIndexing((threadNum - 1)*size, s.length() - K, s, K);
    executor.execute(task);
    tasks[threadNum - 1] = task;
    executor.shutdown();
    executor.awaitTermination(1, TimeUnit.HOURS);
    HashMap<String, TIntArrayList> indexTemp = tasks[0].index;
    tasks[0].index = null;
    for(int n = 1; n < threadNum; n++){
      task = tasks[n];
      for(Map.Entry<String, TIntArrayList> entry: task.index.entrySet()){
        TIntArrayList list = indexTemp.get(entry.getKey());
        if(list == null){
          list = entry.getValue();
          indexTemp.put(entry.getKey(), list);
        }else{
          list.addAll(entry.getValue());
        }
        entry.setValue(null);
      }
      tasks[n].index = null;
    }
    HashMap<String, int[]> index = new HashMap<>();
    for(Map.Entry<String, TIntArrayList> entry: indexTemp.entrySet()){
      index.put(entry.getKey(), entry.getValue().toArray());
      entry.setValue(null);
    }
    return index;
  }

  static HashMap<String, int[]> buildIndex(String seq, int K){
    HashMap<String, TIntArrayList> indexTemp = new HashMap<>();
    int seqLen = seq.length();
    for(int j = 0; j <= seqLen - K; j++){
      String s = seq.substring(j, j + K);
      TIntArrayList l = indexTemp.get(s);
      if(l == null){
        l = new TIntArrayList();
        indexTemp.put(s, l);
      }
      l.add(j);
    }
    HashMap<String, int[]> res = new HashMap<>();
    for(Map.Entry<String, TIntArrayList> entry: indexTemp.entrySet()){
      res.put(entry.getKey(), entry.getValue().toArray());
      entry.setValue(null);
    }
    return res;
  }

}
