// @brief Commissioning analysis of measured data
///
/// @author Adam Mendl <adam.mendl@cvut.cz>


#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <optional>
#include <tuple>
#include <memory>
#include <cstdlib>
#include <cmath>
#include <cfenv>


#include <snfee/snfee.h>
#include <snfee/io/multifile_data_reader.h>

#include <sncabling/om_id.h>
#include <sncabling/gg_cell_id.h>
#include <sncabling/label.h>

#include <snemo/datamodels/raw_event_data.h>
#include <snemo/datamodels/calo_digitized_hit.h>
#include <snemo/datamodels/tracker_digitized_hit.h>

#include <TH1D.h>
#include <TH2D.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TROOT.h>

struct Tracker {
	static const inline double DriftModelA1 		= 828;	/// @brief Betsy Landells: SuperNEMO Drift Model - ns/cm 
	static const inline double DriftModelB1 		= -0.907;	/// @brief Betsy Landells: SuperNEMO Drift Model 
    /// @brief 
    /// 
    /// @param time drift time in ns
    ///
    /// @returns cm
    static double r(const double time) {
        return std::pow(time/DriftModelA1,1/(1-DriftModelB1));
    }
};

/// @brief 
/// @param first_string 
/// @param second_string 
/// @return 
std::vector<std::vector<TH2D*>> permutations(const std::vector<std::string>& first_string,const std::vector<std::string>& second_string) {
    std::vector<std::vector<TH2D*>> retval;
    size_t i = 0;
    size_t i1,i2,i3,i4;
    for(i1 = 0; i1 < second_string.size(); i1++) {
        std::vector<TH2D*> vector2;
        std::vector<std::string> values2;
        for(i2 = 0; i2 < second_string.size(); i2++){
            if(i2==i1) continue;
            values2.push_back(second_string[i2]);
        }
        for(i3 = 0; i3 < first_string.size(); i3++) {
            for(i4 = 0; i4 < first_string.size(); i4++) {
                auto histogram = new TH2D(
                    (first_string[i3] + first_string[i4] + values2[0] + values2[1]).c_str(),
                    (first_string[i3] + "/(" + first_string[i3] + values2[0] + ")" + " to " + first_string[i4] + "/(" + first_string[i4] + values2[1] + ")").c_str(),
                    2*10*20,0,1,2*10*20,0,1
                );
                histogram->GetXaxis()->SetTitle((first_string[i3] + "/(" + first_string[i3] + values2[0] + ")").c_str());
                histogram->GetYaxis()->SetTitle((first_string[i4] + "/(" + first_string[i4] + values2[1] + ")").c_str());
                vector2.push_back(histogram);
            }
        }
        retval.push_back(vector2);
    }
    return retval;
}
/// @brief  
/// @param histograms 
/// @param values1 
/// @param values2 
void fillPermutations(std::vector<std::vector<TH2D*>>& histograms, const std::vector<long double>& values1, const std::vector<long double>& values2){
    size_t i_1 = 0, i_2 = 0;
    size_t i1,i2,i3,i4;
    for(i1 = 0; i1 < values2.size(); i1++) {
        std::vector<long double> values_2;
        for(i2 = 0; i2 < values2.size(); i2++){
            if(i2==i1) continue;
            values_2.push_back(values2[i2]);
        }
        i_2 = 0;
        for(i3 = 0; i3 < values1.size(); i3++) {
            for(i4 = 0; i4 < values1.size(); i4++) {
                histograms[i_1][i_2]->Fill(
                    values1[i3] / (values1[i3] + values_2[0]),
                    values1[i4] / (values1[i4] + values_2[1])
                );
                i_2++;
            }
        }
        i_1++;
    }
}


bool TimestampOK(const snemo::datamodel::timestamp &timestamp) {
	return timestamp.get_ticks() != snfee::data::INVALID_TICKS && timestamp.get_ticks() != 0;
}

/// @brief 
int main(int argc, char** argv) {
    size_t run = std::stoi(argv[1]);

    snfee::initialize();

    snfee::io::multifile_data_reader::config_type reader_cfg;
	reader_cfg.filenames.push_back(std::string("/sps/nemo/snemo/snemo_data/raw_data/v1/RED/delta-tdc-10us/snemo_run-") + std::to_string(run) + "_red.data.gz");
    std::cout<<"Running with run #"<<run<<std::endl;

	// Instantiate a reader
	snfee::io::multifile_data_reader red_source (reader_cfg);

	// Working RED object
	snemo::datamodel::raw_event_data red;
	
	// RED counter
	std::size_t red_counter = 0;


    std::vector<std::string> first_label {
        "T1",
        "T3",
        "Min"
    }; 

    std::vector<std::string> second_label {
        "T2",
        "T4",
        "Max"
    };

    auto histograms = permutations(first_label,second_label);
    std::cout<<histograms.size()<<std::endl;

    while (red_source.has_record_tag()) {
        try{
            DT_THROW_IF(!red_source.record_tag_is(snemo::datamodel::raw_event_data::SERIAL_TAG),
                std::logic_error, "Unexpected record tag '" << red_source.get_record_tag() << "'!");

            red_source.load(red);

            auto tracker_hits       = red.get_tracker_hits();
            auto red_calo_hits 	    = red.get_calo_hits();


            for(const auto& hit: tracker_hits) {
                const auto times = hit.get_times().front();

                if(
                    TimestampOK(times.get_bottom_cathode_time())    &&
                    TimestampOK(times.get_top_cathode_time())       &&
                    TimestampOK(times.get_anode_time(0))            &&
                    TimestampOK(times.get_anode_time(1))            &&
                    TimestampOK(times.get_anode_time(2))            &&
                    TimestampOK(times.get_anode_time(3))            &&
                    TimestampOK(times.get_anode_time(4))
                ) {
                    long double r0 = static_cast<long double>(times.get_anode_time(0).get_ticks())*12.5;

                    std::vector<long double> first_values { 
                        static_cast<long double>(times.get_anode_time(1).get_ticks())*12.5 - r0,
                        static_cast<long double>(times.get_anode_time(3).get_ticks())*12.5 - r0
                    };
                    std::vector<long double> second_values {
                        static_cast<long double>(times.get_anode_time(2).get_ticks())*12.5 - r0,
                        static_cast<long double>(times.get_anode_time(4).get_ticks())*12.5 - r0,
                        static_cast<long double>(std::max(times.get_top_cathode_time().get_ticks(),times.get_bottom_cathode_time().get_ticks()))*12.5 - r0
                    };

                    fillPermutations(histograms,first_values,second_values);
                }
            }


            
        } catch(...) {

        }
        if(red_counter % 10000 == 0)
             std::cout<< "Event #"<<red_counter<<std::endl;
        red_counter++;
    } 

    TFile* file = new TFile((std::to_string(run) + "/multigraph2.root").c_str(),"RECREATE");
    size_t j = 0;
    for(auto&h_ : histograms) {
        auto canvas = new TCanvas();
        size_t q = 1;
        canvas->Divide(2,2);
        for( auto& hist : h_){
            canvas->cd(q);
            hist->Draw("colz");
            q++;
        }

        file->WriteObject(canvas,std::to_string(j).c_str());
        j++;
    }

    delete file;

}
