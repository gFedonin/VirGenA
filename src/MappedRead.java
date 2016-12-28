/**
 * Created by Gennady
 */
class MappedRead implements Comparable{

  int start;
  int end;
  int count;
  Alignment aln;
  byte[] seq;
  String name;
  byte reverse;
  byte[] q;
  byte n;

  MappedRead(){

  }

  MappedRead(int s, int e, int c){
    start = s;
    end = e;
    count = c;
  }

  @Override
  public int compareTo(Object o){
    MappedRead h = (MappedRead) o;
    return h.count - count;
  }

  @Override
  public boolean equals(Object o){
    MappedRead read = (MappedRead) o;
    return read.start == start && read.end == end;
  }

}