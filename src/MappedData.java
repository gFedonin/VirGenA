import java.util.ArrayList;

/**
 * Created by Геннадий on 07.12.2014.
 */
class MappedData {

    ArrayList<MappedRead> mappedReads;
    ArrayList<PairedRead> concordant;
    ArrayList<PairedRead> discordant;
    ArrayList<PairedRead> leftMateMapped;
    ArrayList<PairedRead> rightMateMapped;
    ArrayList<PairedRead> unmapped;

    int[] pairs;
    byte[] isConcordant;
    ArrayList<PairedRead> needToRemap;

    MappedData(){
        mappedReads = new ArrayList<>();
        concordant = new ArrayList<>();
        discordant = new ArrayList<>();
        leftMateMapped = new ArrayList<>();
        rightMateMapped = new ArrayList<>();
        unmapped = new ArrayList<>();
        needToRemap = new ArrayList<>();
    }

}
