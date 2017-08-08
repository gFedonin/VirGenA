/**
 * Created with IntelliJ IDEA.
 * Date: 07.08.15
 */
class Location extends Interval{

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

}
