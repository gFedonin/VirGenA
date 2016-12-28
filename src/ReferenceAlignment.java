import gnu.trove.iterator.TShortObjectIterator;
import gnu.trove.list.array.TShortArrayList;
import gnu.trove.map.hash.TShortObjectHashMap;
import org.jdom2.Document;
import org.jdom2.Element;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileReader;
import java.util.*;

/**
 * Created by Геннадий on 22.03.2015.
 */
class ReferenceAlignment extends Constants{

  ArrayList<Reference> refAlns;
  HashMap<String, AlnPos[]> refAlnIndex;
  HashMap<String, Reference> refSeqs;
  int length;
  int refNum;

  private static ReferenceAlignment instance;

  private ReferenceAlignment(Document document){
    Element element = document.getRootElement().getChild("ReferenceSelector");
    String refAlnPath = element.getChildText("ReferenceMSA");
    String randomModelPath = element.getChild("MapperToMSA").getChildText("RandomModelPath");
    try {
      DataInputStream inputStream = new DataInputStream(new FileInputStream(randomModelPath));
      int k = inputStream.readInt();
      inputStream.close();
      readRefAlns(refAlnPath);
      buildAlnIndex(k);
    }catch (Exception e){
      e.printStackTrace();
    }
  }

  private ReferenceAlignment(int k, String refAlnPath){
    try {
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

  static ReferenceAlignment getInstance(int k, String refAlnPath){
    if(instance == null){
      instance = new ReferenceAlignment(k, refAlnPath);
    }
    return instance;
  }

  private void readRefAlns(String refAlnPath){
    try{
      refAlns = new ArrayList<>();
      refSeqs = new HashMap<>();
      BufferedReader reader = new BufferedReader(new FileReader(refAlnPath));
      String line;
      String name = reader.readLine().substring(1);
      StringBuilder builder = new StringBuilder();
      int refID = 0;
      while((line = reader.readLine()) != null){
        if(line.startsWith(">")){
          Reference seq = new Reference(name, builder.toString());
          seq.refID = refID;
          refID++;
          refAlns.add(seq);
          refSeqs.put(name, seq);
          builder = new StringBuilder();
          name = line.substring(1);
        }else{
          builder.append(line);
        }
      }
      Reference seq = new Reference(name, builder.toString());
      seq.refID = refID;
      refAlns.add(seq);
      refSeqs.put(name, refAlns.get(refAlns.size() - 1));
    }catch(Exception e){
      e.printStackTrace();
    }
    length = refAlns.get(0).aln.length();
    refNum = refAlns.size();
  }

  private class AlnIndexBuilder implements Runnable{

    private Reference seq;
    public short refID;
    public int K;
    public HashMap<String, TShortArrayList> index;

    public AlnIndexBuilder(short refID, Reference seq, int K){
      this.seq = seq;
      this.refID = refID;
      this.K = K;
    }

    @Override
    public void run(){
      index = new HashMap<>();
      String str = seq.seq;
      int[] strToAln = seq.seqToAln;
      for(short i = 0; i <= str.length() - K; i++){
        String s = str.substring(i, i + K);
        if(s.indexOf(GAP) != -1){
          continue;
        }
        TShortArrayList l = index.get(s);
        if(l == null){
          l = new TShortArrayList();
          index.put(s, l);
        }
        l.add((short)strToAln[i]);
      }
    }

  }

  private void buildAlnIndex(int K) throws InterruptedException{
    int threadNum = Runtime.getRuntime().availableProcessors();
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
    HashMap<String, TShortObjectHashMap<TShortArrayList>> indexTemp = new HashMap<>();
    boolean allComplete = false;
    while(!allComplete){
      allComplete = true;
      for(int j = 0; j < threadNum && j < tasks.length; j++){
        if(!threads[j].isAlive()){
          AlnIndexBuilder task = tasks[taskIDs[j]];
          if(task.index == null){
            continue;
          }
          for(Map.Entry<String, TShortArrayList> entry: task.index.entrySet()){
            TShortObjectHashMap<TShortArrayList> map = indexTemp.get(entry.getKey());
            if(map == null){
              map = new TShortObjectHashMap<>();
              indexTemp.put(entry.getKey(), map);
            }
            TShortArrayList list = entry.getValue();
            for(int k = 0; k < list.size(); k++){
              short pos = list.getQuick(k);
              TShortArrayList refList = map.get(pos);
              if(refList == null){
                refList = new TShortArrayList();
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
      Thread.currentThread().sleep(1);
    }
    refAlnIndex = new HashMap<>();
    for(Map.Entry<String, TShortObjectHashMap<TShortArrayList>> entry: indexTemp.entrySet()){
      TShortObjectHashMap<TShortArrayList> v = entry.getValue();
      TShortObjectIterator<TShortArrayList> iter = v.iterator();
      AlnPos[] list = new AlnPos[v.size()];
      for(int j = 0; j < list.length; j++){
        iter.advance();
        TShortArrayList l = iter.value();
        list[j] = new AlnPos();
        list[j].pos = iter.key();
        list[j].refIDs = l.toArray();
      }
      Arrays.sort(list);
      refAlnIndex.put(entry.getKey(), list);
      entry.setValue(null);
    }
  }

}
