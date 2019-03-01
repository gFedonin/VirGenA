import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TIntObjectHashMap;
import org.jdom2.Document;
import org.jdom2.Element;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by Геннадий on 22.03.2015.
 */
class ReferenceAlignment extends Constants{

  ArrayList<Reference> refAlns;
  HashMap<String, int[]> refAlnIndex;
  HashMap<String, Reference> refSeqs;
  int length;
  int refNum;
  int threadNum;

  private static ReferenceAlignment instance;

  private ReferenceAlignment(Document document){
    Element element = document.getRootElement().getChild("ReferenceSelector");
    boolean enabled = Boolean.parseBoolean(element.getChildText("Enabled"));
    if(enabled){
      String refAlnPath = element.getChildText("ReferenceMSA");
      int k = Integer.parseInt(element.getChild("MapperToMSA").getChildText("K"));
      threadNum = Integer.parseInt(document.getRootElement().getChildText("ThreadNumber"));
      if(threadNum == -1){
        threadNum = Runtime.getRuntime().availableProcessors();
      }
      readRefAlns(refAlnPath);
      try{
        buildAlnIndex(k);
      }catch(InterruptedException e){
        e.printStackTrace();
      }
    }
  }

  private ReferenceAlignment(int k, String refAlnPath, int threadNum){
    try {
      this.threadNum = threadNum;
      readRefAlns(refAlnPath);
      buildAlnIndex(k);
    }catch (Exception e){
      e.printStackTrace();
    }
  }

  static ReferenceAlignment getInstance(Document document){
    if(instance == null){
      instance = new ReferenceAlignment(document);
    }
    return instance;
  }

  static ReferenceAlignment getInstance(int k, String refAlnPath, int threadNum){
    if(instance == null){
      instance = new ReferenceAlignment(k, refAlnPath, threadNum);
    }
    return instance;
  }

  private void readRefAlns(String refAlnPath){
    try{
      refAlns = new ArrayList<>();
      refSeqs = new HashMap<>();
      BufferedReader reader = new BufferedReader(new FileReader(refAlnPath));
      String line;
      String name = reader.readLine().substring(1).replace(' ', '_');
      StringBuilder builder = new StringBuilder();
      int refID = 0;
      while((line = reader.readLine()) != null){
        if(line.startsWith(">")){
          Reference ref = new Reference(name, builder.toString());
          ref.refID = refID;
          refID++;
          refAlns.add(ref);
          refSeqs.put(name, ref);
          builder = new StringBuilder();
          name = line.substring(1).replace(' ', '_');
        }else{
          builder.append(line);
        }
      }
      Reference ref = new Reference(name, builder.toString());
      ref.refID = refID;
      refAlns.add(ref);
      refSeqs.put(name, ref);
    }catch(Exception e){
      e.printStackTrace();
    }
    length = refAlns.get(0).aln.length();
    refNum = refAlns.size();
  }

  private class AlnIndexBuilder implements Runnable{

    private Reference reference;
    public int refID;
    public int K;
    public HashMap<String, TIntArrayList> index;

    AlnIndexBuilder(int refID, Reference reference, int K){
      this.reference = reference;
      this.refID = refID;
      this.K = K;
    }

    @Override
    public void run(){
      index = new HashMap<>();
      String str = reference.seq;
      int[] strToAln = reference.seqToAln;
      for(int i = 0; i <= str.length() - K; i++){
        String s = str.substring(i, i + K);
        if(s.indexOf(GAP) != -1){
          continue;
        }
        TIntArrayList l = index.get(s);
        if(l == null){
          l = new TIntArrayList();
          index.put(s, l);
        }
        l.add(strToAln[i]);
      }
    }

  }

  private void buildAlnIndex(int K) throws InterruptedException{
    AlnIndexBuilder[] tasks = new AlnIndexBuilder[refAlns.size()];
    for(short i = 0; i < refAlns.size(); i++){
      Reference reference = refAlns.get(i);
      AlnIndexBuilder task = new AlnIndexBuilder(i, reference, K);
      tasks[i] = task;
    }
    Thread[] threads = new Thread[threadNum];
    int[] taskIDs = new int[threadNum];
    int i = threadNum;
    for(int j = 0; j < threadNum && j < tasks.length; j++){
      Thread thread = new Thread(tasks[j]);
      threads[j] = thread;
      taskIDs[j] = j;
      thread.start();
    }
    HashMap<String, TIntObjectHashMap<TIntArrayList>> indexTemp = new HashMap<>();
    boolean allComplete = false;
    while(!allComplete){
      allComplete = true;
      for(int j = 0; j < threadNum && j < tasks.length; j++){
        if(!threads[j].isAlive()){
          AlnIndexBuilder task = tasks[taskIDs[j]];
          if(task.index == null){
            continue;
          }
          for(Map.Entry<String, TIntArrayList> entry: task.index.entrySet()){
            TIntObjectHashMap<TIntArrayList> map = indexTemp.get(entry.getKey());
            if(map == null){
              map = new TIntObjectHashMap<>();
              indexTemp.put(entry.getKey(), map);
            }
            TIntArrayList list = entry.getValue();
            for(int k = 0; k < list.size(); k++){
              int pos = list.getQuick(k);
              TIntArrayList refList = map.get(pos);
              if(refList == null){
                refList = new TIntArrayList();
                map.put(pos, refList);
              }
              refList.add(task.refID);
            }
            entry.setValue(null);
          }
          task.index = null;
          if(i < tasks.length){
            threads[j] = new Thread(tasks[i]);
            taskIDs[j] = i;
            threads[j].start();
            i++;
            allComplete = false;
          }
        }else{
          allComplete = false;
        }
      }
      Thread.sleep(1);
    }
    refAlnIndex = new HashMap<>();
    for(Map.Entry<String, TIntObjectHashMap<TIntArrayList>> entry: indexTemp.entrySet()){
      TIntObjectHashMap<TIntArrayList> v = entry.getValue();
      TIntObjectIterator<TIntArrayList> iter = v.iterator();
      int[] list = new int[v.size()];
      for(int j = 0; j < list.length; j++){
        iter.advance();
        list[j] = iter.key();
      }
      Arrays.sort(list);
      refAlnIndex.put(entry.getKey(), list);
      entry.setValue(null);
    }
  }

}
