/// @brief Study of tracker timestamps
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

/// @brief Static structure holding information about drift model
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



/// @brief Calculates distance from orthogonal projection
///
/// @param line
/// @param data
/// 
/// @returns
double projection(const std::vector<double> &line,const std::vector<double>& data) {
    size_t i;
    double mod1 = 0., mod2 = 0., dot_product = 0., result = 0.;

    for(double d : line)                mod1 += d*d;
    for(double d : data)                mod2 += d*d;

    mod1 = sqrt(mod1);
    mod2 = sqrt(mod2);

    for(i = 0; i < line.size(); i++)    dot_product += line[i]*data[i];

    dot_product = dot_product / (mod1*mod1);

    for(i = 0; i < line.size(); i++)    result += (data[i] - dot_product*line[i])*(data[i] - dot_product*line[i]);

    return sqrt(result);
}

/// @brief check if value is between a1 and a2
///
///
bool InBounds(long double value, long double a1, long double a2) {
    if(value >= a1 && value <= a2) return true;
    return false;
}


/// @brief Checks if timestamp is OK 
bool TimestampOK(const snemo::datamodel::timestamp &timestamp) {
	return timestamp.get_ticks() != snfee::data::INVALID_TICKS && timestamp.get_ticks() != 0;
}

/// @brief 
int main() {

    // std::vector<double> line = {1.,1.,};
    // std::vector<double> data = {30000.,-30000.};

    // std::cout<<projection(line,data)<<std::endl;
    // std::cin.get();

    size_t run = 807; /// XXX run ID

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

    //std::array<TH1D*,11*142*2> histogramsFirstSignal;
    //std::arayt<TH1D*,11*142*2> histogramsSecondSignal;

    auto r0Histogram            = new TH1D("","",2*10*40,-2*200*12.5,+2*19*200*12.5);
    auto r1Histogram            = new TH1D("","",2*10*400,-2*200*12.5,+2*19*200*12.5);
    auto r2Histogram            = new TH1D("","",2*10*400,-2*200*12.5,+2*19*200*12.5);
    auto r3Histogram            = new TH1D("","",2*10*400,-2*200*12.5,+2*19*200*12.5);
    auto r4Histogram            = new TH1D("","",2*10*400,-2*200*12.5,+2*19*200*12.5);
    auto r3r1PlusHistogram      = new TH1D("","",2*10*400,-2*200*12.5,+2*19*200*12.5);
    auto r4r2PlusHistogram      = new TH1D("","",2*10*400,-2*200*12.5,+2*19*200*12.5);
    auto r5r6PlusHistogram      = new TH1D("","",2*10*400,-2*200*12.5,+2*19*200*12.5);
    auto r1r2PlusHistogram      = new TH1D("","",2*10*400,-2*200*12.5,+2*19*200*12.5);
    auto r3r4PlusHistogram      = new TH1D("","",2*10*400,-2*200*12.5,+2*19*200*12.5);
    auto r5r6MinusHistogram     = new TH1D("r5r6Minus","r5 - r6",2*10*400,-2*200*12.5,+2*19*200*12.5);
    auto r1r2MinusHistogram     = new TH1D("r1r2Minus","r1 - r2",2*10*400,-2*200*12.5,+2*19*200*12.5);
    auto speedHistogram         = new TH1D("speed","speed",2*10*400,1./90000.,1./30000.);

    
    auto r                      = new TH1D("r","r",40,0,8);

    auto CreateTH2D = [](const char* name, const char * title, const char* title_x, const char * title_y) {
        auto retval = new TH2D(name,title,2*10*20,-2*200*12.5,+2*19*200*12.5,2*10*20,-2*200*12.5,+2*19*200*12.5);

        retval->GetXaxis()->SetTitle(title_x);
        retval->GetYaxis()->SetTitle(title_y);
        
        return retval;
    };

    auto cathodeAnodeCorrelationHistogram           = CreateTH2D("cathodeAnodeCorrelationHistogram", "Correlation between sum of cathode and sum of anode signals (ns)","#frac{R_1 + R_2 + R_3 + R_4}{2}","R_5 + R_6");
    auto firstSignalCorrelationHistogram            = CreateTH2D("firstSignalCorrelationHistogram","Correlation between cathode and anode for first signal (ns)","#frac{R_1 + R_3}{2}","min(R_5,R_6)");
    auto secondSignalCorrelationHistogram           = CreateTH2D("secondSignalCorrelationHistogram","Correlation between cathode and anode for second signal (ns)","#frac{R_2 + R_4}{2}","max(R_5,R_6)");
    auto anodeFirstSecondCorrelationHistogram       = CreateTH2D("anodeFirstSecondCorrelation","Correlation between first and second signal for anode","#frac{R_1 + R_3}{2}","#frac{R_2 + R_4}{2}");
    auto r1r3CorrelationHistogram                   = CreateTH2D("r1r3CorrelationHistogram","Correlation between R1 and R3","R1","R3");
    auto r2r4CorrelationHistogram                   = CreateTH2D("r2r4CorrelationHistogram","Correlation between R2 and R4","R2","R4");

    auto zSpeedCorrelationHistogram                 = new TH2D("zSpeedCorrelationHistogram","Correlation between z and t_5 + t_6", 2*10*20,  -2*200*12.5,+2*19*200*12.5,  2*10*20, -1., 1.);


    

    auto r3r1MinusHistogram     = new TH1D("","",2*10*400,-2*200*12.5,+2*19*200*12.5);
    auto r4r2MinusHistogram     = new TH1D("","",2*10*400,-2*200*12.5,+2*19*200*12.5);
    auto meanToFirstCathode     = new TH1D("","",2*10*400,-2*200*12.5,+2*19*200*12.5);
    auto meanToSecondCathode    = new TH1D("","",2*10*400,-2*200*12.5,+2*19*200*12.5);

    TH1D* histogramFirstSignal      = new TH1D("","Correlation of cathode and anode - first signal",300,-150.*12.5,150.*12.5);
    TH1D* histogramSecondSignal     = new TH1D("","Correlation of cathode and anode - second signal",300,-150.*12.5,150.*12.5);
    TH1D* propagationTimeHistogram  = new TH1D("","",7000,0,7000);
    TH1D* differenceHistogram       = new TH1D("","",20000,0,20000);
    
    size_t calohitsConsidered = 0;
    while (red_source.has_record_tag()) {
        if(red_counter%10000 == 0) std::cout<<red_counter<<"\n"<<std::flush;
        try{
            DT_THROW_IF(!red_source.record_tag_is(snemo::datamodel::raw_event_data::SERIAL_TAG),
                std::logic_error, "Unexpected record tag '" << red_source.get_record_tag() << "'!");

            red_source.load(red);

            auto tracker_hits       = red.get_tracker_hits();
            auto red_calo_hits 	    = red.get_calo_hits();

            
            size_t calohitsConsideredInEvent = 0;
            for(const auto& tracker_hit : tracker_hits) {
                auto times = tracker_hit.get_times().front();

                // check sanity
                if(
                    TimestampOK(times.get_bottom_cathode_time())    &&
                    TimestampOK(times.get_top_cathode_time())       &&
                    TimestampOK(times.get_anode_time(0))            &&
                    TimestampOK(times.get_anode_time(1))            &&
                    TimestampOK(times.get_anode_time(2))            &&
                    TimestampOK(times.get_anode_time(3))            &&
                    TimestampOK(times.get_anode_time(4))
                ) {
                    auto first                  = std::min(times.get_bottom_cathode_time(),times.get_top_cathode_time());
                    auto second                 = std::max(times.get_bottom_cathode_time(),times.get_top_cathode_time());

                    histogramFirstSignal        ->Fill((times.get_anode_time(1).get_ticks() + times.get_anode_time(3).get_ticks() - 2*(first.get_ticks()))*12.5);
                    histogramSecondSignal       ->Fill((times.get_anode_time(2).get_ticks() + times.get_anode_time(4).get_ticks() - 2*(second.get_ticks()))*12.5);
                    propagationTimeHistogram    ->Fill((second.get_ticks() + first.get_ticks() - 2*times.get_anode_time(0).get_ticks())*12.5);                    
                    differenceHistogram         ->Fill((second.get_ticks() - first.get_ticks())*12.5);

                    size_t index                                = 0;
                    size_t index_best                           = -1;
                    size_t smallest_ticks                       = -1;

                    calohitsConsideredInEvent = 0;
                    int good_calohits = 0;
                    for(const auto& hit : red_calo_hits) {
                        if(!hit.is_high_threshold()) {
                            continue;
                        }
                        else
                            good_calohits++;
                        calohitsConsideredInEvent++;
                        if(good_calohits>1) break;


                        long double hit_ns = static_cast<long double>(hit.get_reference_time().get_ticks())*6.25 - 300.;
                        if(
                            index_best == -1 || (
                                hit_ns <= static_cast<long double>(times.get_anode_time(0).get_ticks())*12.5 &&
                                static_cast<long double>(red_calo_hits[index_best].get_reference_time().get_ticks())*6.25- 300. < hit_ns
                            )
                        ) {
                            index_best = index;
                        }
                        if(
                            smallest_ticks == -1 ||
                            hit_ns < static_cast<long double>(red_calo_hits[smallest_ticks].get_reference_time().get_ticks())*6.25 - 300.
                        ) {
                            smallest_ticks = index;
                        }
                        index++;

                        
                        
                    }
                    if(index_best == -1 && smallest_ticks == -1 ) continue;
                    long double s0 = index_best != -1? static_cast<long double>(red_calo_hits[index_best].get_reference_time().get_ticks())*6.25 - 300. : static_cast<long double>(red_calo_hits[smallest_ticks].get_reference_time().get_ticks())*6.25 - 300.;
                    if(good_calohits>1) continue;

                    long double R0                                                  = static_cast<long double>(times.get_anode_time(0).get_ticks())*12.5-s0;
                    long double R1                                                  = static_cast<long double>(times.get_anode_time(1).get_ticks())*12.5-s0;
                    long double R2                                                  = static_cast<long double>(times.get_anode_time(2).get_ticks())*12.5-s0;
                    long double R3                                                  = static_cast<long double>(times.get_anode_time(3).get_ticks())*12.5-s0;
                    long double R4                                                  = static_cast<long double>(times.get_anode_time(4).get_ticks())*12.5-s0;

                    long double R5                                                  = first .get_ticks()*12.5 - s0;
                    long double R6                                                  = second.get_ticks()*12.5 - s0;

                    long double R5_real                                             = times.get_top_cathode_time().get_ticks()*12.5 - s0;
                    long double R6_real                                             = times.get_bottom_cathode_time().get_ticks()*12.5 - s0;


                    std::vector<double> line = {1.,1.};
                    std::vector<double> data = {(R1 + R3)/2.,R5};


                    double zRelative                                                = (R5_real - R6_real)/(R5_real + R6_real - 2*R0);

                    // if(!InBounds(1/(R5_real + R6_real - 2*R0),17E-6,20E-6))
                    // // if(!InBounds(1/(R5_real + R6_real - 2*R0),21E-6,24E-6))
                    //      continue;

                    // if(!(abs(zRelative) < 0.7 && abs(zRelative) > 0.3))
                    // if(!(abs(zRelative) < 0.6 && abs(zRelative) > 0.4))
                        // continue;
                    // if(projection(line,data) > 2000.)
                    //     continue;

    // R2 < 24500. ||
                    // if( Tracker::r(R0) > 2.2)
                    //     continue;


                    zSpeedCorrelationHistogram                                      ->Fill(R5_real + R6_real - 2*R0,(R5_real - R6_real)/(R5_real + R6_real - 2*R0));
                    speedHistogram                                                  ->Fill(1/(R5_real + R6_real - 2*R0));


                    r->Fill(Tracker::r(R0));
                    // std::cout<<Tracker::r(R0)<<std::endl;
                    // std::cin.get();


                    r0Histogram                                                     ->Fill(R0);
                    r1Histogram                                                     ->Fill(R1);
                    r2Histogram                                                     ->Fill(R2);
                    r3Histogram                                                     ->Fill(R3);
                    r4Histogram                                                     ->Fill(R4);

                    cathodeAnodeCorrelationHistogram                                ->Fill((R1 + R2 + R3 + R4)/2,R5 + R6);
                    firstSignalCorrelationHistogram                                 ->Fill((R1 + R3)/2,R5);
                    secondSignalCorrelationHistogram                                ->Fill((R2+R4)/2,R6);
                    anodeFirstSecondCorrelationHistogram                            ->Fill((R1 + R3)/2,(R2 + R4)/2);
                    r1r3CorrelationHistogram                                        ->Fill(R1,R3);
                    r2r4CorrelationHistogram                                        ->Fill(R2,R4);

                    r5r6MinusHistogram                                              ->Fill((times.get_bottom_cathode_time().get_ticks() - times.get_top_cathode_time().get_ticks())*12.5);
                    r1r2MinusHistogram                                              ->Fill(R2-R1);

                    // std::cout<<times.get_anode_time(3).get_ticks()*12.5<<" "<<times.get_anode_time(1).get_ticks()*12.5<<std::endl<<s0<<std::endl;
                    // std::cin.get();

                    r3r1PlusHistogram           ->Fill(static_cast<long double>(times.get_anode_time(3).get_ticks() + times.get_anode_time(1).get_ticks())*12.5 - 2*s0);
                    r4r2PlusHistogram           ->Fill(static_cast<long double>(times.get_anode_time(4).get_ticks() + times.get_anode_time(2).get_ticks())*12.5 - 2*s0);

                    r3r1MinusHistogram          ->Fill((times.get_anode_time(3) - times.get_anode_time(1)).get_ticks()*12.5);
                    r4r2MinusHistogram          ->Fill((times.get_anode_time(4) - times.get_anode_time(2)).get_ticks()*12.5);

                    meanToFirstCathode          ->Fill((static_cast<long double>(times.get_anode_time(3).get_ticks() + times.get_anode_time(1).get_ticks())*12.5)/2 - first.    get_ticks()*12.5);
                    meanToSecondCathode         ->Fill((static_cast<long double>(times.get_anode_time(4).get_ticks() + times.get_anode_time(2).get_ticks())*12.5)/2 - second.   get_ticks()*12.5);

                    r5r6PlusHistogram           ->Fill(static_cast<long double>(times.get_top_cathode_time().get_ticks() + times.get_bottom_cathode_time().get_ticks())*12.5 - 2*s0);

                    r1r2PlusHistogram           ->Fill((times.get_anode_time(1).get_ticks() + times.get_anode_time(2).get_ticks())*12.5 - 2*s0);
                    r3r4PlusHistogram           ->Fill((times.get_anode_time(3).get_ticks() + times.get_anode_time(4).get_ticks())*12.5 - 2*s0);

                }
            }

            calohitsConsidered+=calohitsConsideredInEvent;
        } catch (...) {
            std::cerr<<"Error occured in event #"<<red_counter<<'\t'<<std::flush;
        }
        red_counter++;
    }
    std::cout<<"Processed "<<red_counter<< " events"<<std::endl;

    TFile* file = new TFile((std::to_string(run) + "/output.root").c_str(),"RECREATE"); // XXX save root file

    auto save = [&](TH1* h,std::string title,std::string name) {
        auto c = new TCanvas();
        h->SetNameTitle(name.c_str(),title.c_str());
        h->Draw();
        c->SaveAs((std::to_string(run) + '/' + (name) + ".pdf").c_str());
        file->WriteObject(h,name.c_str());

    };

    save(r0Histogram,"0","0");
    save(r1Histogram,"1","1");
    save(r2Histogram,"2","2");
    save(r3Histogram,"3","3");
    save(r4Histogram,"4","4");

    save(r1r2PlusHistogram,"R1 + R2 (nano s)","1plus2");
    save(r3r4PlusHistogram,"R3 + R4 (nano s)","3plus4");

    save(r3r1MinusHistogram,"R3 - R1 (nano s)","3minus1");
    save(r4r2MinusHistogram,"R4 - R2 (nano s)","4minus2");
    save(r3r1PlusHistogram,"R3 + R1 (nano s)","3plus1");
    save(r4r2PlusHistogram,"R4 + R2 (nano s)","4plus2");
    save(r5r6PlusHistogram,"R5 + R6 (nano s)","5plus6");

    file->WriteObject(cathodeAnodeCorrelationHistogram,cathodeAnodeCorrelationHistogram->GetName());
    file->WriteObject(firstSignalCorrelationHistogram,firstSignalCorrelationHistogram->GetName());
    file->WriteObject(secondSignalCorrelationHistogram,secondSignalCorrelationHistogram->GetName());
    file->WriteObject(anodeFirstSecondCorrelationHistogram,anodeFirstSecondCorrelationHistogram->GetName());
    file->WriteObject(r1r3CorrelationHistogram,r1r3CorrelationHistogram->GetName());
    file->WriteObject(r2r4CorrelationHistogram,r2r4CorrelationHistogram->GetName());
    file->WriteObject(r,r->GetName());
    file->WriteObject(r1r2MinusHistogram,r1r2MinusHistogram->GetName());
    file->WriteObject(r5r6MinusHistogram,r5r6MinusHistogram->GetName());
    file->WriteObject(zSpeedCorrelationHistogram,zSpeedCorrelationHistogram->GetName());
    file->WriteObject(speedHistogram,speedHistogram->GetName());

    // save(histogramFirstSignal,"canvas1","canvas1");
    // save(histogramSecondSignal,"canvas2","canvas1");

    save(meanToFirstCathode,"mean to first cathode (nano s)","MeanToFirstCathode");
    save(meanToSecondCathode,"mean to second cathode (nano s)","MeanToSecondCathode");


    delete file;
    std::cout<<"Calohits considered "<<calohitsConsidered<<std::endl;

    return 0;

}
