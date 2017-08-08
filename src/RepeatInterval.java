/**
 * Created by Gennady on 26.05.2017.
 */
public class RepeatInterval extends Interval implements Comparable{

  int assemblyID;
  int startDest;
  int endDest;

  RepeatInterval(int s, int e, int id, int sD, int eD){
    start = s;
    end = e;
    assemblyID = id;
    startDest = sD;
    endDest = eD;
  }

  @Override
  public int compareTo(Object o){
    RepeatInterval interval = (RepeatInterval) o;
    return start - interval.start;
  }

  @Override
  public boolean equals(Object o){
    RepeatInterval interval = (RepeatInterval) o;
    return start == interval.start && end == interval.end &&
        assemblyID == interval.assemblyID && startDest == interval.startDest &&
        endDest == interval.endDest;
  }
}
