import org.jdom2.Document;
import org.jdom2.Element;
import org.jdom2.JDOMException;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;

/**
 * Created by Gennady on 20.01.2016.
 */
class Postprocessor extends Constants{

  private float minIdentity;
  private float minCoverage;
  private int minLen;
  private String outPath;
  private boolean debug;
  private int threadNum;

  Postprocessor(Document document){
    Element element = document.getRootElement().getChild("Postprocessor");
    debug = Boolean.parseBoolean(element.getChildText("Debug"));
    outPath = document.getRootElement().getChildText("OutPath");
    minLen = Integer.parseInt(element.getChildText("MinFragmentLength"));
    minIdentity = Float.parseFloat(element.getChildText("MinIdentity"));
    threadNum = Integer.parseInt(document.getRootElement().getChildText("ThreadNumber"));
    minCoverage = Float.parseFloat(element.getChildText("MinFragmentCoverage"));
  }

  private void makeBlastDB(String[] assemblies, Reference[] refSeqs) throws IOException, InterruptedException{
    for(int i = 0; i < assemblies.length; i++){
      String seq = assemblies[i];
      String name = refSeqs[i].name;
      Process p = Runtime.getRuntime().exec("makeblastdb -in - -parse_seqids -dbtype nucl -title " +
        name + " -out " + outPath + name);
      BufferedWriter streamWriter =
          new BufferedWriter(new OutputStreamWriter(p.getOutputStream()));
      streamWriter.write(">" + name + "\n");
      streamWriter.write(seq);
      streamWriter.close();
      p.waitFor();
    }
    Process p = Runtime.getRuntime().exec("makeblastdb -in - -parse_seqids -dbtype nucl " +
        "-title refSeqDB -out " + outPath + "refSeqDB");
    BufferedWriter streamWriter =
        new BufferedWriter(new OutputStreamWriter(p.getOutputStream()));
    for(Reference seq: refSeqs){
      streamWriter.write(">" + seq.name + "\n");
      streamWriter.write(seq.seq);
      streamWriter.newLine();
    }
    streamWriter.close();
    p.waitFor();
  }

  private void deleteBlastDB(Reference[] refSeqs){
    for(Reference seq: refSeqs){
      File file = new File(outPath + seq.name + ".nhr");
      file.delete();
      file = new File(outPath + seq.name + ".nin");
      file.delete();
      file = new File(outPath + seq.name + ".nog");
      file.delete();
      file = new File(outPath + seq.name + ".nsd");
      file.delete();
      file = new File(outPath + seq.name + ".nsi");
      file.delete();
      file = new File(outPath + seq.name + ".nsq");
      file.delete();
    }
    File file = new File(outPath + "refSeqDB.nhr");
    file.delete();
    file = new File(outPath + "refSeqDB.nin");
    file.delete();
    file = new File(outPath + "refSeqDB.nog");
    file.delete();
    file = new File(outPath + "refSeqDB.nsd");
    file.delete();
    file = new File(outPath + "refSeqDB.nsi");
    file.delete();
    file = new File(outPath + "refSeqDB.nsq");
    file.delete();
  }

  private float refIdentity(HSP interval, Reference seq1, Reference seq2, int[] seqToAln){
    int identity = 0;
    int len = 0;
    int start = seqToAln[interval.startS];
    int end;
    if(interval.endS == seqToAln.length){
      end = seq1.aln.length();
    }else{
      end = seqToAln[interval.endS];
    }
    for(int i = start; i < end; i++){
      char c1 = seq1.aln.charAt(i);
      char c2 = seq2.aln.charAt(i);
      if(c1 == GAP){
        if(c2 != GAP){
          len ++;
        }
      }else{
        if(c2 == c1){
          identity ++;
        }
        len ++;
      }
    }
    if(debug){
      System.out.println(seq1.name + " vs " + seq2.name);
      System.out.println(interval.startS + " " + interval.endS + " id = " + (float) identity/len);
    }
    return (float)identity/len;
  }

  private int chooseRef(HSP interval, String assembly, Reference ref1, Reference ref2, int readNum1, int readNum2) throws IOException, JDOMException, InterruptedException{
    String hitSeq = assembly.substring(interval.startQ, interval.endQ);
//    Process p = Runtime.getRuntime().exec(
//        "blastn -outfmt 5 -db refSeqDB", new String[]{"BLASTDB=" + outPath}, null);
    Process p = Runtime.getRuntime().exec(
            "blastn -num_threads " + threadNum + " -outfmt 10 -db refSeqDB", new String[]{"BLASTDB=" + outPath}, null);
    BufferedWriter streamWriter =
        new BufferedWriter(new OutputStreamWriter(p.getOutputStream()));
    streamWriter.write(hitSeq);
    streamWriter.close();
    ArrayList<BlastHit> list;
    if(debug){
      BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
      BufferedWriter writer = new BufferedWriter(new FileWriter(outPath + ref1.name +
          "_" + interval.startQ + "_" + interval.endQ + ".csv"));
      String line;
      while((line = reader.readLine()) != null){
        writer.write(line);
        writer.newLine();
      }
      writer.close();
      reader = new BufferedReader(new InputStreamReader(p.getErrorStream()));
      while((line = reader.readLine()) != null){
        System.out.println(line);
      }
      p.waitFor();
      list = BlastParser.parseBlastCSV(new FileInputStream(outPath + ref1.name +
          "_" + interval.startQ + "_" + interval.endQ + ".csv"));
    }else{
      list = BlastParser.parseBlastCSV(p.getInputStream());
      BufferedReader reader = new BufferedReader(new InputStreamReader(p.getErrorStream()));
      String line;
      while((line = reader.readLine()) != null){
        System.out.println(line);
      }
      p.waitFor();
    }
    //ArrayList<BlastHit> list = parseBlastXML(p.getInputStream());
//    ArrayList<BlastHit> list = BlastParser.parseBlastCSV(p.getInputStream());
//    p.waitFor();
    HSP bestInt1 = null;
    HSP bestInt2 = null;
    for(BlastHit hit: list){
      if(hit.hitID.equals(ref1.name)){
        int maxLen = 0;
        for(HSP interval1: hit.hsps){
          if(interval1.endQ - interval1.startQ > maxLen){
            maxLen = interval1.endQ - interval1.startQ;
            bestInt1 = interval1;
          }
        }
      }else if(hit.hitID.equals(ref2.name)){
        int maxLen = 0;
        for(HSP interval2: hit.hsps){
          if(interval2.endQ - interval2.startQ > maxLen){
            maxLen = interval2.endQ - interval2.startQ;
            bestInt2 = interval2;
          }
        }
      }
    }
    if(bestInt1 == null){
      if(bestInt2 != null){
        return 2;
      }
    }else{
      if(bestInt2 == null){
        return 1;
      }else{
        if(readNum1 <= readNum2){
          if(refIdentity(bestInt2, ref1, ref2, ref2.seqToAln) < minIdentity){
            return 2;
          }else{
            return 12;
          }
        }else{
          if(refIdentity(bestInt1, ref1, ref2, ref1.seqToAln) < minIdentity){
            return 1;
          }else{
            return 12;
          }
        }
      }
    }
    return 12;
  }

  private ArrayList<HSP> getSimilarFragments(String seqQ, int[] contigEndsQ, String refNameQ, String refNameS,
                                             int[] contigEndsS) throws IOException, InterruptedException, JDOMException{
    ArrayList<HSP> hsps = new ArrayList<>();
//    Process p = Runtime.getRuntime().exec(
//        "blastn -outfmt 5 -db " + dbName, new String[]{"BLASTDB=" + outPath}, null);
    Process p = Runtime.getRuntime().exec(
            "blastn -num_threads " + threadNum + " -outfmt 10 -db " + refNameS, new String[]{"BLASTDB=" + outPath}, null);
    BufferedWriter streamWriter =
        new BufferedWriter(new OutputStreamWriter(p.getOutputStream()));
    streamWriter.write(seqQ);
    streamWriter.close();
    ArrayList<BlastHit> list;
    if(debug){
      BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
      BufferedWriter writer = new BufferedWriter(new FileWriter(outPath + refNameQ + "-" + refNameS + ".csv"));
      String line;
      while((line = reader.readLine()) != null){
        writer.write(line);
        writer.newLine();
      }
      writer.close();
      reader = new BufferedReader(new InputStreamReader(p.getErrorStream()));
      while((line = reader.readLine()) != null){
        System.out.println(line);
      }
      p.waitFor();
      list = BlastParser.parseBlastCSV(new FileInputStream(outPath + refNameQ + "-" + refNameS + ".csv"));
    }else{
      list = BlastParser.parseBlastCSV(p.getInputStream());
      p.waitFor();
    }
    //ArrayList<BlastHit> list = parseBlastXML(p.getInputStream());
    //ArrayList<BlastHit> list = BlastParser.parseBlastCSV(p.getInputStream());
    //p.waitFor();
    if(list.isEmpty()){
      return hsps;
    }
    for(BlastHit blastHit: list){
      for(HSP hsp: blastHit.hsps){
        if(hsp.identity >= minIdentity){
          if(Math.min(hsp.endQ - hsp.startQ, hsp.endS - hsp.startS) >= minLen){
            hsps.add(hsp);
          }else{
            int startOfContigQ = getContigStart(hsp.startQ, contigEndsQ);
            int endOfContigQ = getContigEnd(hsp.endQ, contigEndsQ);
            int startOfContigS = getContigStart(hsp.startS, contigEndsS);
            int endOfContigS = getContigEnd(hsp.endS, contigEndsS);
            float covQ = (float)(hsp.endQ - hsp.startQ)/(endOfContigQ - startOfContigQ);
            float covS = (float)(hsp.endS - hsp.startS)/(endOfContigS - startOfContigS);
            if(covQ >= minCoverage || covS >= minCoverage){
              hsps.add(hsp);
            }
          }
        }
      }
    }
    return hsps;
  }

  private static int getContigStart(int pos, int[] contigEnds){
    int startOfContig = 0;
    for(int contigEnd : contigEnds){
      if(pos < contigEnd){
        return startOfContig;
      }else{
        startOfContig = contigEnd;
      }
    }
    return startOfContig;
  }

  private static int getContigEnd(int pos, int[] contigEnds){
    int endOfContig = 0;
    for(int contigEnd : contigEnds){
      if(pos <= contigEnd){
        return contigEnd;
      }
    }
    return endOfContig;
  }


  private ArrayList<RepeatInterval>[] getRawIntervals(String[] assemblies, int[][] contigEnds, Reference[] selectedRefs, int[] readNums) throws InterruptedException, JDOMException, IOException{
    ArrayList<RepeatInterval>[] intervals = new ArrayList[assemblies.length];
    for(int i = 0; i < intervals.length; i++){
      intervals[i] = new ArrayList<>();
    }
    for(int i = 0; i < assemblies.length; i++){
      String seq_i = assemblies[i];
      int[] contigEnds_i = contigEnds[i];
      Reference ref_i = selectedRefs[i];
      for(int j = i + 1; j < selectedRefs.length; j++){
        Reference ref_j = selectedRefs[j];
        ArrayList<HSP> similarFragments = getSimilarFragments(seq_i, contigEnds_i, ref_i.name,
            ref_j.name, contigEnds[j]);
        for(HSP interval: similarFragments){
          int refID = chooseRef(interval, seq_i, ref_i, ref_j, readNums[i], readNums[j]);
          switch(refID){
            case 1:
              intervals[j].add(new RepeatInterval(interval.startS, interval.endS, i, interval.startQ, interval.endQ));
              break;
            case 2:
              intervals[i].add(new RepeatInterval(interval.startQ, interval.endQ, j, interval.startS, interval.endS));
              break;
          }
        }
      }
    }
    return intervals;
  }

  private ArrayList<RepeatInterval>[] filterIntervals(ArrayList<RepeatInterval>[] intervals, int[] readNums){
    ArrayList<RepeatInterval>[] res = new ArrayList[intervals.length];
    for(int j = 0; j < intervals.length; j++){
      ArrayList<RepeatInterval> list = intervals[j];
      if(list.isEmpty()){
        res[j] = list;
        continue;
      }
      TreeSet<RepeatInterval> sortedSet = new TreeSet<>(list);
      ArrayList<RepeatInterval> fIntervals = new ArrayList<>();
      RepeatInterval pivot = sortedSet.pollFirst();
      while(!sortedSet.isEmpty()){
        RepeatInterval curr = sortedSet.pollFirst();
        if(pivot.end > curr.start){
          // intersection
          if(pivot.assemblyID == curr.assemblyID){
            if(curr.end > pivot.end){
              pivot.end = curr.end;
            }
          }else{
            if(pivot.end == curr.end && pivot.start == curr.start){
              // this is the same fragment similar to multiple assemblies
              if(readNums[pivot.assemblyID] < readNums[curr.assemblyID]){
                pivot = curr;
              }
            }else{
              if(pivot.end - pivot.start >= curr.end - curr.start){
                // cut or merge curr to pivot
                if(curr.end > pivot.end){
                  curr.start = pivot.end;
                  //if(curr.end - curr.start >= minLen){
                  //this fragment should be cut anyway
                  sortedSet.add(curr);
                  //}
                }
              }else{
                // cut or merge pivot to curr
                pivot.end = curr.start;
                //if(pivot.end - pivot.start >= minLen){
                //this fragment should be cut anyway
                fIntervals.add(pivot);
                //}
                pivot = curr;
              }
            }
          }
        }else{
          fIntervals.add(pivot);
          pivot = curr;
        }
      }
      fIntervals.add(pivot);
      res[j] = fIntervals;
    }
    return res;
  }

  private void setDest(int assemblyID, int num, RepeatInterval[][] intervalsArr, ArrayList<RepeatInterval>[] res){
    RepeatInterval interval = intervalsArr[assemblyID][num];
    RepeatInterval[] destList = intervalsArr[interval.assemblyID];
    RepeatInterval temp = new RepeatInterval(interval.startDest, interval.endDest, 0, 0, 0);
    int i = Arrays.binarySearch(destList, temp);
    if(i >= 0){
      //found interval starting from the same coordinate
      RepeatInterval dest = destList[i];
      res[assemblyID].add(interval);
      if(dest.end < interval.endDest){
        RepeatInterval right = new RepeatInterval(interval.start + dest.end - dest.start,
            interval.end, interval.assemblyID, dest.end, interval.endDest);
        res[assemblyID].add(right);
        interval.end = right.start;
        interval.endDest = right.startDest;
      }
      interval.assemblyID = dest.assemblyID;
    }else{
      int insertPos = -i - 1;
      if(insertPos == destList.length){
        // might be intersection with prev
        if(destList.length > 0){
          RepeatInterval prev = destList[destList.length - 1];
          res[assemblyID].add(interval);
          if(interval.startDest < prev.end){
            // intersects with prev
            if(interval.endDest > prev.end){
              // longer than prev
              RepeatInterval right = new RepeatInterval(interval.start + prev.end - interval.startDest,
                  interval.end, interval.assemblyID, prev.end, interval.endDest);
              res[assemblyID].add(right);
              interval.end = right.start;
              interval.endDest = right.startDest;
            }
            interval.assemblyID = prev.assemblyID;
          }
        }
      }else{
        if(insertPos == 0){
          // might be intersection with next
          RepeatInterval next = destList[insertPos];
          if(interval.endDest > next.end){
            RepeatInterval left = new RepeatInterval(interval.start,
                interval.end + next.start - interval.endDest, interval.assemblyID,
                interval.startDest, next.start);
            res[assemblyID].add(left);
            interval.startDest = left.endDest;
            interval.start = left.end;
            res[assemblyID].add(interval);
            RepeatInterval right = new RepeatInterval(interval.end - interval.endDest + next.end,
                interval.end, interval.assemblyID, next.end, interval.endDest);
            res[assemblyID].add(right);
            interval.end = right.start;
            interval.endDest = right.startDest;
            interval.assemblyID = next.assemblyID;
          }else{
            if(interval.endDest > next.start){
              RepeatInterval left = new RepeatInterval(interval.start,
                  interval.end + next.start - interval.endDest, interval.assemblyID,
                  interval.startDest, next.start);
              res[assemblyID].add(left);
              interval.startDest = left.endDest;
              interval.start = left.end;
              interval.assemblyID = next.assemblyID;
            }
            res[assemblyID].add(interval);
          }
        }else{
          // might be intersection with both prev and next
          RepeatInterval prev = destList[insertPos - 1];
          if(interval.startDest < prev.end){
            // intersects with prev
            if(interval.endDest > prev.end){
              // longer than prev
              RepeatInterval left = new RepeatInterval(interval.start,
                  interval.end - interval.endDest + prev.end, prev.assemblyID,
                  interval.startDest, prev.end);
              res[assemblyID].add(left);
              interval.start = left.end;
              interval.startDest = left.endDest;
            }else{
              // inside the prev
              interval.assemblyID = prev.assemblyID;
              return;
            }
          }
          res[assemblyID].add(interval);
          RepeatInterval next = destList[insertPos];
          if(interval.endDest > next.start){
            //intersects next
            RepeatInterval right = new RepeatInterval(interval.end - interval.endDest + next.start,
                interval.end, next.assemblyID, next.start, interval.endDest);
            res[assemblyID].add(right);
            interval.end = right.start;
            interval.endDest = right.startDest;
          }
        }
      }
    }
  }

  private ArrayList<RepeatInterval>[] updateIntervals(ArrayList<RepeatInterval>[] intervals){
    ArrayList<RepeatInterval>[] res = null;
    outer:
    while(true){
      res = new ArrayList[intervals.length];
      RepeatInterval[][] intervalsArr = new RepeatInterval[intervals.length][];
      for(int i = 0; i < intervals.length; i++){
        res[i] = new ArrayList<>();
        intervalsArr[i] = intervals[i].toArray(new RepeatInterval[intervals[i].size()]);
      }
      for(int i = 0; i < intervals.length; i++){
        for(int j = 0; j < intervalsArr[i].length; j++){
          setDest(i, j, intervalsArr, res);
        }
      }
      for(int i = 0; i < intervals.length; i++){
        if(intervals[i].size() != res[i].size()){
          intervals = res;
          continue outer;
        }
      }
      for(int i = 0; i < intervals.length; i++){
        if(!intervals[i].equals(res[i])){
          intervals = res;
          continue outer;
        }
      }
      break;
    }
    return res;
  }

  ArrayList<RepeatInterval>[] findRegionsToCut(String[] assemblies, int[][] contigEnds, Reference[] selectedRefs, int[] readNums) throws IOException, InterruptedException, JDOMException{
    makeBlastDB(assemblies, selectedRefs);
    ArrayList<RepeatInterval>[] intervals = getRawIntervals(assemblies, contigEnds, selectedRefs, readNums);
    ArrayList<RepeatInterval>[] res = filterIntervals(intervals, readNums);
    res = updateIntervals(res);
    deleteBlastDB(selectedRefs);
    if(debug){
      BufferedWriter writer = new BufferedWriter(new FileWriter(outPath + "intervals.txt"));
      for(int i = 0; i < res.length; i++){
        writer.write(selectedRefs[i].name + "\n");
        for(RepeatInterval interval: res[i]){
          writer.write(Integer.toString(interval.start) + " " +
              Integer.toString(interval.end) + " " + selectedRefs[interval.assemblyID].name + "\n");
        }
      }
      writer.close();
    }
    return res;
  }

}
