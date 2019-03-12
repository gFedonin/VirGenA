/**
 * Created by Gennady
 */
class Alignment{

  /**
   * Alignment score
   */
  int score;

  /**
   * Aligned sequence #1
   */
  byte[] sequence1;

  /**
   * Alignment start location in sequence #1
   */
  int start1;
  int end1;

  /**
   * Alignment start location in sequence #2
   */
  int start2;
  int end2;

  /**
   * Count of identical bases
   */
  int identity;

  int length;

}
