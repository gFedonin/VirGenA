/**
 * Created with IntelliJ IDEA.
 * Date: 07.08.15
 */
class Location extends Interval implements Comparable{

  int startIndex;
  int endIndex;
  int count;

  Location(){

  }

  Location(Location st){
    start = st.start;
    startIndex = st.startIndex;
    end = st.end;
    endIndex = st.endIndex;
    count = st.count;
  }

  @Override
  public int compareTo(Object o){
    Location l = (Location) o;
    if(start == l.start){
      if(end == l.end){
        return 0;
      }
      return (end < l.end)?-1:1;
    }
    return (start < l.start)?-1:1;
  }
}
