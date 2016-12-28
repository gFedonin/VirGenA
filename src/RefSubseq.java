/**
 * Created with IntelliJ IDEA.
 * Date: 22.07.15
 */
class RefSubseq implements Comparable{

  int refID;
  int start;
  int end;

  @Override
  public int compareTo(Object o){
    RefSubseq ref = (RefSubseq) o;
    return (ref.end - ref.start) - (end - start);
  }

  RefSubseq(int refID, ReferenceAlignment refAlignment, int s, int e){
    Reference ref = refAlignment.refAlns.get(refID);
    this.refID = refID;
    start = ref.alnToSeqStart(s);
    if(start == -1){
      end = -1;
      return;
    }
    end = ref.alnToSeqEnd(e);
    if(end == -1){
      start = -1;
    }
  }

  public RefSubseq(int r, int s, int e){
    refID = r;
    start = s;
    end = e;
  }


}
