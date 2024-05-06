/*
  val3dity 

  Copyright (c) 2011-2023, 3D geoinformation research group, TU Delft

  This file is part of val3dity.

  val3dity is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  val3dity is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with val3dity.  If not, see <http://www.gnu.org/licenses/>.

  For any information or further details about the use of val3dity, contact
  Hugo Ledoux
  <h.ledoux@tudelft.nl>
  Faculty of the Built Environment & Architecture
  Delft University of Technology
  Julianalaan 134, Delft 2628BL, the Netherlands
*/

#include "input.h"
#include "Primitive.h"
#include "Surface.h"
#include "MultiSurface.h"
#include "CompositeSurface.h"
#include "Solid.h"
#include "Feature.h"

#include "GenericObject.h"

#include <tclap/CmdLine.h>
#include <time.h>  
#include "nlohmann-json/json.hpp"
#include <boost/filesystem.hpp>

using namespace std;
using namespace val3dity;
using json = nlohmann::json;

std::string VAL3DITY_VERSION = "2.4.1b0";


std::string print_summary_validation(std::vector<Feature*>& lsFeatures, IOErrors& ioerrs);
std::string unit_test(std::vector<Feature*>& lsFeatures, IOErrors& ioerrs);
void        write_report_json(json& jr, std::string report);



class MyOutput : public TCLAP::StdOutput
{
public:
  
  virtual void usage(TCLAP::CmdLineInterface& c)
  {
    std::cout << "USAGE:" << std::endl;
    std::cout << "\tval3dity [OPTION] myinput.json" << std::endl;
    std::cout << "ALLOWED FORMATS:" << std::endl;
    std::cout << "\tCityJSON" << std::endl;
    std::cout << "\tCityJSONSeq" << std::endl;
    std::cout << "\ttu3djson" << std::endl;
    std::cout << "\tJSON-FG" << std::endl;
    std::cout << "\tOBJ" << std::endl;
    std::cout << "\tOFF" << std::endl;
    std::cout << "\tIndoorGML" << std::endl;
    std::cout << "OPTIONS:" << std::endl;
    std::list<TCLAP::Arg*> args = c.getArgList();
    for (TCLAP::ArgListIterator it = args.begin(); it != args.end(); it++) {
      if ( ((*it)->getName() == "ishell") || ((*it)->getName() == "output_off") )
        continue;
      if ((*it)->getFlag() == "")
        std::cout << "\t--" << (*it)->getName() << std::endl;
      else
        std::cout << "\t-" << (*it)->getFlag() << ", --" << (*it)->getName() << std::endl;
      std::cout << "\t\t" << (*it)->getDescription() << std::endl;
    }
    std::cout << "SOME EXAMPLES:" << std::endl;
    
    std::cout << "\tval3dity input.city.json" << std::endl;
    std::cout << "\t\tValidate each City Object and each 3D primitive in input.json (CityJSON file)" << std::endl;
    std::cout << "\t\tand print a summary" << std::endl;
    
    std::cout << "\tval3dity input.json --report /home/elvis/temp/myreport.json" << std::endl;
    std::cout << "\t\tValidate each 3D primitive in input.json (tu3djson file)" << std::endl;
    std::cout << "\t\tand output a detailed JSON report '/home/elvis/temp/myreport.json';" << std::endl;
    std::cout << "\t\tbrowse that report at http://geovalidation.bk.tudelft.nl/val3dity/browser/" << std::endl;
    
    std::cout << "\tval3dity input.city.json --verbose" << std::endl;
    std::cout << "\t\tAll details of the validation are printed out" << std::endl;
    
    std::cout << "\tval3dity input.obj" << std::endl;
    std::cout << "\t\tValidate the geometries in input.obj as if they were an ISO19107 Solid (default)" << std::endl;
    
    std::cout << "\tval3dity input.off -p MultiSurface" << std::endl;
    std::cout << "\t\tValidate the geometries in input.off as an ISO19107 MultiSurface" << std::endl;
    
    std::cout << "\tval3dity myindoorgml.gml --snap_tol 0.1" << std::endl;
    std::cout << "\t\tThe vertices in myindoorgml.gml closer than 0.1unit are snapped together" << std::endl;
  }

  virtual void failure(TCLAP::CmdLineInterface& c, TCLAP::ArgException& e)
  {
    std::cout << "ERROR: " << e.argId() << std::endl; 
    std::cout << "       " << e.error() << endl;
    std::cout << std::endl;
    std::cout << "USAGE: val3dity [OPTION] myinput.json" << std::endl;
    std::cout << std::endl;
    std::cout << "Try `val3dity --help` for more options" << std::endl;
    exit(1);
  }

};

class LicensePrint : public TCLAP::Visitor
{
  public:
  
    LicensePrint() : Visitor() {};
    void visit() 
    {
      std::string thelicense =
      "\nval3dity\n\n"
      "Copyright (c) 2011-2023  3D geoinformation research group, TU Delft\n\n"
      "val3dity is free software: you can redistribute it and/or modify\n"
      "it under the terms of the GNU General Public License as published by\n"
      "the Free Software Foundation, either version 3 of the License, or\n"
      "(at your option) any later version.\n\n"
      "val3dity is distributed in the hope that it will be useful,\n"
      "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
      "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
      "GNU General Public License for more details.\n\n"
      "A copy of the GNU General Public License is available at\n"
      "<http://www.gnu.org/licenses/> or\n"
      "<https://github.com/tudelft3d/val3dity/blob/master/LICENSE\n\n"
      "For any information or further details about the use of val3dity, contact:\n"
      "Hugo Ledoux \n"
      "<h.ledoux@tudelft.nl>\n"
      "Faculty of the Built Environment & Architecture\n"
      "Delft University of Technology\n"
      "Julianalaan 134, Delft 2628BL, the Netherlands\n";
      std::cout << thelicense << std::endl;  
      exit(0); 
    };
};


class ListErrorsPrint : public TCLAP::Visitor
{
  public:
  
    ListErrorsPrint() : Visitor() {};
    void visit() 
    {
      for (auto& e: ALL_ERRORS) {
        std::cout << e.first << " -- " << e.second << std::endl;
      }
      std::cout << std::endl;
      std::string url = "https://val3dity.readthedocs.io/en/" + VAL3DITY_VERSION + "/errors/";
      std::cout << "Explanations and examples of each error at " << url << std::endl;
      exit(0); 
    };
};

int main(int argc, char* const argv[])
{
  IOErrors ioerrs;
  std::streambuf* savedBufferCLOG;
  std::ofstream mylog;

  //-- tclap options
  std::vector<std::string> primitivestovalidate;
  primitivestovalidate.push_back("Solid");  
  primitivestovalidate.push_back("CompositeSurface");   
  primitivestovalidate.push_back("MultiSurface");   
  TCLAP::ValuesConstraint<std::string> primVals(primitivestovalidate);

  TCLAP::CmdLine cmd("Allowed options", ' ', VAL3DITY_VERSION);
  MyOutput my;
  cmd.setOutput(&my);
  try {
    TCLAP::UnlabeledValueArg<std::string>   inputfile(
                                              "inputfile", 
                                              "allowed formats: CityJSON, CityJSONSeq, tu3djson, JSON-FG, IndoorGML, OBJ, or OFF",
                                              true, 
                                              "", 
                                              "string");
    TCLAP::MultiArg<std::string>            ishellfiles(
                                              "", 
                                              "ishell", 
                                              "one interior shell (in POLY format only; more than one possible)",
                                              false, 
                                              "string");
    TCLAP::ValueArg<std::string>            report("r",
                                              "report",
                                              "output report in JSON format",
                                              false,
                                              "",
                                              "string");    
    TCLAP::ValueArg<std::string>            primitives("p",
                                              "primitive",
                                              "geometric type of the input: <Solid|CompositeSurface|MultiSurface> (default=Solid)",
                                              false,
                                              "Solid",
                                              &primVals);
    TCLAP::SwitchArg                        verbose("",
                                              "verbose",
                                              "verbose output",
                                              false);
    TCLAP::SwitchArg                        license("l",
                                              "license",
                                              "see the software license",
                                              false,
                                              new LicensePrint() );
    TCLAP::SwitchArg                        listerrors("",
                                              "listerrors",
                                              "list all the possible errors",
                                              false,
                                              new ListErrorsPrint() );
    TCLAP::SwitchArg                        unittests("",
                                              "unittests",
                                              "unit tests output",
                                              false);
    TCLAP::SwitchArg                        ignore204("",
                                              "ignore204",
                                              "ignore error 204",
                                              false);    
    TCLAP::ValueArg<std::string>            output_off("",
                                              "output_off",
                                              "output each shell/surface in OFF format",
                                              false,
                                              "",
                                              "string");        
    TCLAP::ValueArg<double>                 snap_tol("",
                                              "snap_tol",
                                              "tolerance for snapping vertices (default=0.001)",
                                              false,
                                              0.001,
                                              "double");
    TCLAP::ValueArg<double>                 overlap_tol("",
                                              "overlap_tol",
                                              "tolerance for testing overlap CompositeSolids and BuildingParts (default=0.0)",
                                              false,
                                              -1,
                                              "double");
    TCLAP::ValueArg<double>                 planarity_d2p_tol("",
                                              "planarity_d2p_tol",
                                              "tolerance for planarity distance_to_plane (default=0.01)",
                                              false,
                                              0.01,
                                              "double");
    TCLAP::ValueArg<double>                 planarity_n_tol("",
                                              "planarity_n_tol",
                                              "tolerance for planarity based on normals deviation (default=20.0degree)",
                                              false,
                                              20.0,
                                              "double");

    cmd.add(planarity_d2p_tol);
    cmd.add(planarity_n_tol);
    cmd.add(snap_tol);
    cmd.add(overlap_tol);
    cmd.add(verbose);
    cmd.add(primitives);
    cmd.add(ignore204);
    cmd.add(unittests);
    cmd.add(output_off);
    cmd.add(inputfile);
    cmd.add(listerrors);
    cmd.add(license);
    cmd.add(ishellfiles);
    cmd.add(report);
    cmd.parse( argc, argv );

    std::string licensewarning =
    "---\nval3dity Copyright (c) 2011-2023, 3D geoinformation research group, TU Delft  \n"
    "This program comes with ABSOLUTELY NO WARRANTY.\n"
    "This is free software, and you are welcome to redistribute it\n"
    "under certain conditions; for details run val3dity with the '--license' option.\n---";
    std::cout << licensewarning << std::endl;

    //-- vector with Features: CityObject, GenericObject, 
    //-- or IndoorModel (or others in the future)
    std::vector<Feature*> lsFeatures;
    
    InputTypes inputtype = OTHER;
    std::string extension = inputfile.getValue().substr(inputfile.getValue().find_last_of(".") + 1);
    if ( (extension == "gml") || (extension == "GML") || (extension == "xml") || (extension == "XML") ) {
      inputtype = GML;
      ioerrs.set_input_file_type("GML");
    }
    else if ( (extension == "poly") || (extension == "POLY") ) {
      inputtype = POLY;    
      ioerrs.set_input_file_type("POLY");
    }
    else if ( (extension == "json") || (extension == "JSON") ) {
      inputtype = JSON;
      ioerrs.set_input_file_type("JSON");
    }
    else if ( (extension == "jsonl") || (extension == "JSONL") ) {
      inputtype = JSONL;
      ioerrs.set_input_file_type("JSONL");
    }
    else if ( (extension == "obj") || (extension == "OBJ") ) {
      inputtype = OBJ;
      ioerrs.set_input_file_type("OBJ");
    }
    else if ( (extension == "off") || (extension == "OFF") ) {
      inputtype = OFF;
      ioerrs.set_input_file_type("OFF");
    }

    //-- if verbose == false then log to a file
    if (verbose.getValue() == false)
    {
      savedBufferCLOG = clog.rdbuf();
      mylog.open("val3dity.log");
      std::clog.rdbuf(mylog.rdbuf());
    }

    //-- no negative snap_tol value
    if (snap_tol.getValue() < 0) 
    {
      ioerrs.add_error(903, "snap_tol cannot be negative");
    }

    
    if (inputtype == OTHER) {
      std::stringstream ss;
      ss << "Format of file " << inputfile.getValue() << " not supported (based on its extension).";
      ioerrs.add_error(904, ss.str());
      ioerrs.set_input_file_type("UNKNOWN");
    }

    Primitive3D prim3d = ALL;
    if ( (inputtype == JSON) || (inputtype == GML) )
    {
      prim3d = ALL;
    }
    else {
      if (primitives.getValue() == "Solid")
        prim3d = SOLID;
      else if (primitives.getValue() == "MultiSurface")
        prim3d = MULTISURFACE;
      else if (primitives.getValue() == "CompositeSurface")
        prim3d = COMPOSITESURFACE;
      if ((prim3d == COMPOSITESURFACE) && (ishellfiles.getValue().size() > 0))
        ioerrs.add_error(903, "POLY files having inner shells cannot be validated as CompositeSurface (only Solids)");
    }

    if (ioerrs.has_errors() == false)
    {
      if (inputtype == GML)
      {
        try
        {
          read_file_gml(inputfile.getValue(), 
                        lsFeatures,
                        ioerrs, 
                        snap_tol.getValue());
          if (ioerrs.has_errors() == true) {
            std::cout << "Errors while reading the input file, aborting." << std::endl;
            std::cout << ioerrs.get_report_text() << std::endl;
          }
          if (ishellfiles.getValue().size() > 0)
          {
            std::cout << "No inner shells allowed when GML file used as input." << std::endl;
            ioerrs.add_error(901, "No inner shells allowed when GML file used as input.");
          }
        }
        catch (int e)
        {
          if (e == 901)
            ioerrs.add_error(901, "Invalid GML structure, or that particular construction of GML is not supported. Please report at https://github.com/tudelft3d/val3dity/issues and provide the file.");
        }
      }
      else if (inputtype == JSON)
      {
        read_file_json(inputfile.getValue(), 
                       lsFeatures,
                       ioerrs, 
                       snap_tol.getValue());
        if (ioerrs.has_errors() == true) {
          std::cout << "Errors while reading the input file, aborting." << std::endl;
          std::cout << ioerrs.get_report_text() << std::endl;
        }
        if (ishellfiles.getValue().size() > 0)
        {
          std::cout << "No inner shells allowed when JSON file used as input." << std::endl;
          ioerrs.add_error(901, "No inner shells allowed when JSON file used as input.");
        }
      }
      else if (inputtype == JSONL)
      {
        read_file_jsonl(inputfile.getValue(), 
                        lsFeatures,
                        ioerrs, 
                        snap_tol.getValue());
        if (ioerrs.has_errors() == true) {
          std::cout << "Errors while reading the input file, aborting." << std::endl;
          std::cout << ioerrs.get_report_text() << std::endl;
        }
        if (ishellfiles.getValue().size() > 0)
        {
          std::cout << "No inner shells allowed when JSONL file used as input." << std::endl;
          ioerrs.add_error(901, "No inner shells allowed when JSONL file used as input.");
        }
      }      
      else if (inputtype == POLY)
      {
        std::cout << "Reading file: " << inputfile.getValue() << std::endl;
        std::ifstream infile(inputfile.getValue().c_str(), std::ifstream::in);
        if (!infile) {
          ioerrs.add_error(901, "Input file not found.");
        } else {
          GenericObject* o = new GenericObject("none");
          Surface* sh = parse_poly(infile, 0, ioerrs);
          if ( (ioerrs.has_errors() == false) & (prim3d == SOLID) )
          {
            Solid* s = new Solid;
            s->set_oshell(sh);
            int sid = 1;
            for (auto ifile : ishellfiles.getValue())
            {
              std::ifstream if2(ifile.c_str(), std::ifstream::in);
              Surface* sh = parse_poly(if2, sid, ioerrs);
              if (ioerrs.has_errors() == false)
              {
                s->add_ishell(sh);
                sid++;
              }
            }
            if (ioerrs.has_errors() == false)
              o->add_primitive(s);
          }
          else if ( (ioerrs.has_errors() == false) & (prim3d == COMPOSITESURFACE) )
          {
            CompositeSurface* cs = new CompositeSurface;
            cs->set_surface(sh);
            if (ioerrs.has_errors() == false)
              o->add_primitive(cs);
          }
          else if ( (ioerrs.has_errors() == false) & (prim3d == MULTISURFACE) )
          {
            MultiSurface* ms = new MultiSurface;
            ms->set_surface(sh);
            if (ioerrs.has_errors() == false)
              o->add_primitive(ms);
          }
          lsFeatures.push_back(o);      
        }
      }
      else if (inputtype == OFF)
      {
        std::cout << "Reading file: " << inputfile.getValue() << std::endl;
        std::ifstream infile(inputfile.getValue().c_str(), std::ifstream::in);
        if (!infile) {
          ioerrs.add_error(901, "Input file not found.");
        } else {
          GenericObject* o = new GenericObject("none");
          Surface* sh = parse_off(infile, 0, ioerrs, snap_tol.getValue());
          if ( (ioerrs.has_errors() == false) & (prim3d == SOLID) )
          {
            Solid* s = new Solid;
            s->set_oshell(sh);
            if (ioerrs.has_errors() == false)
              o->add_primitive(s);
          }
          else if ( (ioerrs.has_errors() == false) & (prim3d == COMPOSITESURFACE) )
          {
            CompositeSurface* cs = new CompositeSurface;
            cs->set_surface(sh);
            if (ioerrs.has_errors() == false)
              o->add_primitive(cs);
          }
          else if ( (ioerrs.has_errors() == false) & (prim3d == MULTISURFACE) )
          {
            MultiSurface* ms = new MultiSurface;
            ms->set_surface(sh);
            if (ioerrs.has_errors() == false)
              o->add_primitive(ms);
          }
          lsFeatures.push_back(o);
        }
      }    
      else if (inputtype == OBJ)
      {
        std::cout << "Reading file: " << inputfile.getValue() << std::endl;
        std::ifstream infile(inputfile.getValue().c_str(), std::ifstream::in);
        if (!infile) {
          ioerrs.add_error(901, "Input file not found.");
        } else {
          parse_obj(infile, lsFeatures, prim3d, ioerrs, snap_tol.getValue());
          if (ioerrs.has_errors() == true) {
            std::cout << "Errors while reading the input file, aborting." << std::endl;
            std::cout << ioerrs.get_report_text() << std::endl;
          }
          if (ishellfiles.getValue().size() > 0)
          {
            std::cout << "No inner shells allowed when GML file used as input." << std::endl;
            ioerrs.add_error(901, "No inner shells allowed when GML file used as input.");
          }
        }
      }
    }

    double planarity_n_tol_updated = planarity_n_tol.getValue();
    if (ioerrs.has_errors() == false) 
    {
      //-- report on parameters used
      if (ignore204.getValue() == true)
        planarity_n_tol_updated = 180.0;
      std::cout << "Parameters used for validation:" << std::endl;
      if (snap_tol.getValue() < 0)
        std::cout << "   snap_tol"    << setw(22)  << "0.001" << std::endl;
      else
        std::cout << "   snap_tol"    << setw(22)  << snap_tol.getValue() << std::endl;
      std::cout << "   planarity_d2p_tol"     << setw(13)  << planarity_d2p_tol.getValue() << std::endl;
      std::cout << "   planarity_n_tol"       << setw(15) << planarity_n_tol_updated << std::endl;
      if (overlap_tol.getValue() < 1e-8)
        std::cout << "   overlap_tol" << setw(19)  << "none" << std::endl;
      else
        std::cout << "   overlap_tol" << setw(19)  << overlap_tol.getValue() << std::endl;
      std::cout << std::endl;
    }
    
    //-- now the validation starts
    if ( (lsFeatures.empty() == false) && (ioerrs.has_errors() == false) )
    {
      int i = 1;
      std::cout << "Validation of " << lsFeatures.size() << " feature(s):" << std::endl;
      for (auto& f : lsFeatures)
      {
        if ( (i % 10 == 0) && (verbose.getValue() == false) )
          printProgressBar(100 * (i / double(lsFeatures.size())));
        i++;
        f->validate(planarity_d2p_tol.getValue(), planarity_n_tol_updated, overlap_tol.getValue());
      }
      if (verbose.getValue() == false)
        printProgressBar(100);
    }

    //-- if error 901 then ignore what was read, it can't be validated
    //-- and is confusing for users to see a valid/invalid while nothing was done...
    if (ioerrs.has_specific_error(901) == true) {
      // std::cout << "ERROR 901" << std::endl;
      lsFeatures.clear();
    }

    //-- summary of the validation
    std::cout << "\n" << print_summary_validation(lsFeatures, ioerrs) << std::endl;        

    //-- output shells/surfaces in OFF format
    if (output_off.getValue() != "") 
    {
      std::cout << std::endl << std::endl;
      boost::filesystem::path outpath(output_off.getValue());
      if (boost::filesystem::exists(outpath.parent_path()) == false)
        std::cout << "Error OFF output: file " << outpath << " impossible to create, wrong path." << std::endl;
      else {
        if (boost::filesystem::exists(outpath) == false)
          boost::filesystem::create_directory(outpath);
        std::cout << "OFF files saved to: " << outpath << std::endl;
        for (auto& f : lsFeatures) {
          int noprim = 0;
          for (auto& p : f->get_primitives()) {
            std::string theid = f->get_id() + "." + std::to_string(noprim);
            if (p->get_type() == SOLID)
            {
              boost::filesystem::path outfile = outpath / (theid + ".0.off");
              std::ofstream o(outfile.string());
              Solid* ts = dynamic_cast<Solid*>(p);
              o << ts->get_off_representation(0) << std::endl;                                
              o.close();
              for (int i = 1; i <= ts->num_ishells(); i++)
              {
                outfile = outpath / (theid + "." + std::to_string(i) + ".off");
                o.open(outfile.string());
                o << ts->get_off_representation(1) << std::endl;                                
                o.close();
              }
            }
            else if (p->get_type() == MULTISURFACE)
            {
              boost::filesystem::path outfile = outpath / (theid + ".off");
              std::ofstream o(outfile.string());
              MultiSurface* ts = dynamic_cast<MultiSurface*>(p);
              o << ts->get_off_representation() << std::endl;                                
              o.close();
            }
            else if (p->get_type() == COMPOSITESURFACE)
            {
              boost::filesystem::path outfile = outpath / (theid + ".off");
              std::ofstream o(outfile.string());
              CompositeSurface* ts = dynamic_cast<CompositeSurface*>(p);
              o << ts->get_off_representation() << std::endl;                                
              o.close();
            }            
            else {
              std::cout << "OFF OUTPUT: these primitive types are not supported (yet). Sorry." << std::endl;
            }
            noprim++;
          }
        }
      }
      std::cout << std::endl;
    }

    //-- output report in JSON 
    if (report.getValue() != "") 
    {
      //-- save the json report in memory first
      json jr = get_report_json(inputfile.getValue(),
                                lsFeatures,
                                VAL3DITY_VERSION,
                                snap_tol.getValue(),
                                overlap_tol.getValue(),
                                planarity_d2p_tol.getValue(),
                                planarity_n_tol_updated,
                                ioerrs);
      if (report.getValue() != "")
        write_report_json(jr, report.getValue());
    }
    else
      std::cout << "==> The validation report wasn't saved, use option '--report'." << std::endl;

    //-- unittests 
    if (unittests.getValue() == true)
      std::cout << "\n" << unit_test(lsFeatures, ioerrs) << std::endl;

    if (verbose.getValue() == false)
    {
      clog.rdbuf(savedBufferCLOG);
      mylog.close();
    }
    return(0);
  }
  catch (TCLAP::ArgException &e) 
  {
    std::cout << "ERROR: " << e.error() << " for arg " << e.argId() << std::endl;
    return(0);
  }
}


void write_report_json(json& jr, std::string report)
{
  boost::filesystem::path outpath(report);
  if (boost::filesystem::exists(outpath.parent_path()) == false) {
    boost::filesystem::path fullpath(boost::filesystem::current_path() / report);
    if (boost::filesystem::exists(fullpath.parent_path()) == false) {
      std::cout << "Error: file " << outpath << " impossible to create, wrong path." << std::endl;
      return;
    }
    else
      outpath = fullpath;
  }
  if (outpath.extension() != ".json")
    outpath += ".json";
  std::ofstream o(outpath.string());
  o << jr.dump(2) << std::endl;                                
  std::cout << "Validation report saved to " << boost::filesystem::canonical(outpath) << std::endl;
  std::cout << "Browse its content:" << std::endl;
  std::cout << "==>http://geovalidation.bk.tudelft.nl/val3dity/browser/\n" << std::endl;
}


std::string unit_test(std::vector<Feature*>& lsFeatures, IOErrors& ioerrs)
{
  std::stringstream ss;
  ss << std::endl;
  std::set<int> theerrors;
  if (ioerrs.has_errors() == true)
  {
    for (auto& each : ioerrs.get_unique_error_codes())
      theerrors.insert(each);
  }
  for (auto& f : lsFeatures) {
    for (auto& code : f->get_unique_error_codes())
        theerrors.insert(code);
    for (auto& p : f->get_primitives())
      for (auto& code : p->get_unique_error_codes())
        theerrors.insert(code);
  }
  if (theerrors.size() > 0)
  {
    ss << "@INVALID ";
    for (auto each : theerrors)
      ss << each << " ";
  }
  else {
    ss << "@VALID";
  }
  ss << std::endl;
  return ss.str();
}


std::string print_summary_validation(std::vector<Feature*>& lsFeatures, IOErrors& ioerrs)
{
  std::stringstream ss;
  ss << std::endl;
  int noprim = 0;
  for (auto& o : lsFeatures)
    for (auto& p : o->get_primitives())
      noprim++;
  //-- overview of errors
  std::map<int,int> errors_f; //-- features
  for (auto& f : lsFeatures)
    for (auto& code : f->get_unique_error_codes())
      if (code > 600)
        errors_f[code] = 0;
  for (auto& f : lsFeatures)
    for (auto& code : f->get_unique_error_codes())
      if (code > 600)
        errors_f[code] += 1;
  std::map<int,int> errors_p; //-- primitives
  for (auto& f : lsFeatures)
    for (auto& p : f->get_primitives())
      for (auto& code : p->get_unique_error_codes())
        errors_p[code] = 0;
  for (auto& f : lsFeatures)
    for (auto& p : f->get_primitives())
      for (auto& code : p->get_unique_error_codes())
        errors_p[code] += 1;

  ss << "+++++++++++++++++++ SUMMARY +++++++++++++++++++" << std::endl;
  if ( (errors_f.size() > 0) || (errors_p.size() > 0) || (ioerrs.has_errors() == true) )
    ss << "INVALID :(" << std::endl;
  else
    ss << "VALID :)" << std::endl;
  ss << "+++++" << std::endl;
  ss << "Input file type:" << std::endl;
  std::string ft = ioerrs.get_input_file_type();
  ss << "  " << ft << std::endl;
  ss << "+++++" << std::endl;
  int fInvalid = 0;
  for (auto& f : lsFeatures)
  {
    if (f->is_valid() == false)
      fInvalid++;
  }
  ss << "Total # of Features: " << setw(10) << lsFeatures.size() << std::endl;
  float percentage;
  if (lsFeatures.size() == 0)
    percentage = 0;
  else
    percentage = 100 * (fInvalid / float(lsFeatures.size()));
  ss << "  # valid: " << setw(20) << lsFeatures.size() - fInvalid;
  if (lsFeatures.size() == 0)
    ss << " (" << 0 << "%)" << std::endl;
  else
    ss << std::fixed << setprecision(1) << " (" << 100 - percentage << "%)" << std::endl;
  ss << "  # invalid: " << setw(18) << fInvalid;
  ss << std::fixed << setprecision(1) << " (" << percentage << "%)" << std::endl;
  std::set<std::string> thetypes;
  for (auto& f : lsFeatures)
    thetypes.insert(f->get_type());
  if (thetypes.empty() == false)
  {
    ss << "Types:" << std::endl;
    for (auto& each : thetypes)
      ss << "  " << each << std::endl;
  }
  ss << "+++++" << std::endl;
  ss << "Total # of primitives: " << setw(8) << noprim << std::endl;
  int bValid = 0;
  for (auto& f : lsFeatures)
    for (auto& p : f->get_primitives())
      if (p->is_valid() == true)
        bValid++;
  if (noprim  == 0)
    percentage = 0;
  else
    percentage = 100 * ((noprim - bValid) / float(noprim));
  ss << "  # valid: " << setw(20) << bValid;
  if (noprim == 0)
    ss << " (" << 0 << "%)" << std::endl;
  else
    ss << std::fixed << setprecision(1) << " (" << 100 - percentage << "%)" << std::endl;
  ss << "  # invalid: " << setw(18) << (noprim - bValid);
  ss << std::fixed << setprecision(1) << " (" << percentage << "%)" << std::endl;
  std::set<int> theprimitives;
  for (auto& f : lsFeatures)
    for (auto& p : f->get_primitives())
      theprimitives.insert(p->get_type());
  if (theprimitives.empty() == false)
  {
    ss << "Types:" << std::endl;
    for (auto& each : theprimitives)
    {
      ss << "  ";
      switch(each)
      {
        case 0: ss << "Solid"             << std::endl; break;
        case 1: ss << "CompositeSolid"    << std::endl; break;
        case 2: ss << "MultiSolid"        << std::endl; break;
        case 3: ss << "CompositeSurface"  << std::endl; break;
        case 4: ss << "MultiSurface"      << std::endl; break;
        case 5: ss << "GeometryTemplate"  << std::endl; break;
        case 9: ss << "ALL"               << std::endl; break;
      }
    }
  }
  if ( (errors_f.size() > 0) || (errors_p.size() > 0) || (ioerrs.has_errors() == true) )
  {
    ss << "+++++" << std::endl;
    ss << "Errors present:" << std::endl;
    //-- primitives
    for (auto e : errors_p)
    {
      ss << "  " << e.first << " -- " << ALL_ERRORS[e.first] << std::endl;
      ss << setw(11) << e.second;
      ss << " primitive(s)";
      ss << std::endl;
    }
    for (auto e : errors_f)
    {
      ss << "  " << e.first << " -- " << ALL_ERRORS[e.first] << std::endl;
      ss << setw(11) << e.second;
      ss << " feature(s)";
      ss << std::endl;
    }
    for (auto& e : ioerrs.get_unique_error_codes())
    {
      ss << "  " << e << " -- " << ALL_ERRORS[e] << std::endl;
    }
  }
  ss << "+++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  return ss.str();
}




