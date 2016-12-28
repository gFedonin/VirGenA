import java.util.ArrayList;

/**
 * Created by Геннадий on 09.12.2014.
 */
abstract class ConsensusBuilder extends Constants{

  public abstract String buildConsensus(Reference genome, ArrayList<PairedRead> reads);

}
