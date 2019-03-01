import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.*;

public class VirGenA{

  public static void main(String[] args){
    ArgumentParser parser = ArgumentParsers.newFor("java -jar VirGenA.jar").build();
    Subparsers subparsers = parser.addSubparsers().title("subcommands").help("description:").dest("command").metavar("COMMAND");

    Subparser parserMapper = subparsers.addParser("map").help("maps given reads to reference");
    Mapper.addParameters(parserMapper);
    Subparser parserRMBuilder = subparsers.addParser("model").help("[deprecated] builds model for given reference or reference MSA");
    RandomModelBuilder.addParameters(parserRMBuilder);
    Subparser parserRefFinder = subparsers.addParser("type").help("chooses the best reference subset from given MSA");
    ReferenceFinder.addParameters(parserRefFinder);
    Subparser parserAssemble = subparsers.addParser("assemble").help("assembles the given paired reads with assist " +
        "of given reference or reference MSA");
    RefBasedAssembler.addParameters(parserAssemble);
    Subparser parserTAB = subparsers.addParser("tab").help("generates TAB file and corresponding LIB file " +
        "needed to run SSPACE scaffolder using VirGenA assembly with BAM file.");
    GenerateTabFile.addParameters(parserTAB);
    if(args.length == 0){
      parser.printUsage();
      return;
    }
    try{
      Namespace parsedArgs = parser.parseArgs(args);
      String command = parsedArgs.getString("command");
      switch(command){
        case "map":
          Mapper.run(parsedArgs);
          break;
        case "assemble":
          RefBasedAssembler.run(parsedArgs);
          break;
        case "type":
          ReferenceFinder.run(parsedArgs);
          break;
        case "model":
          RandomModelBuilder.run(parserRMBuilder, parsedArgs);
          break;
        case "tab":
          GenerateTabFile.run(parsedArgs);
          break;
        default:
          parser.printUsage();
      }
    }catch(ArgumentParserException e){
      parser.handleError(e);
    }
  }

}
