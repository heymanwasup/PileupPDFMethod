#include <vector>
#include <iostream>
#include <string>
#include <map>
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMath.h"

using namespace std;


const double nSlices = 4691; // 4691 * 149.2 ns ~ 700e3 ns
constexpr int nPuEnergyBins = 62;


class PileupHolder {
  public:
    PileupHolder();
    void SetETHist(TH2* th2_et);
    void SetTimeSlice(int n_t_slice);
    std::pair<double,double> GetPDF_d(const double& e_x,const double& e_y);
    std::pair<double,double> GetPDF_t(const double& e_x,const double& e_y,const double& e_z);

    std::pair<double,double> GetPDF_d(const int& e1Bin,const int& e2Bin);
    std::pair<double,double> GetPDF_t(const int& e1Bin,const int& e2Bin,const int& e3Bin);
    int GetNSlices() const {return dist->GetXaxis()->GetNbins();}
    double GetBinCenter(const int& nbin);
    int FindBin(const double& center);
    TH2 * GetHist() {return dist;}

  private:
    TH2 * dist;
    vector<vector<std::pair<float,float>>> data;
    int n_t_slice;
    double DT;

    double eWidth;
    double eStart;
};

double PileupHolder::GetBinCenter(const int& nbin){
    return eStart + eWidth * (nbin - 0.5);
}

int PileupHolder::FindBin(const double& center) {
    return int((center - eStart) / eWidth) + 1;
}

PileupHolder::PileupHolder(){
    DT = 1e-4;
}

void PileupHolder::SetETHist(TH2* th2_et) {
    dist = th2_et;
    eWidth = dist->GetYaxis()->GetBinWidth(1);
    eStart = dist->GetYaxis()->GetBinLowEdge(1);
    
    data.clear();
    int Nx = dist->GetXaxis()->GetNbins();
    int Ny = dist->GetYaxis()->GetNbins();    
    for(int nx=0;nx<Nx;nx++){
        data.push_back(vector<std::pair<float,float>>());
        for(int ny=0;ny<Ny;ny++){
            data[nx].push_back(std::make_pair(dist->GetBinContent(nx+1,ny+1),dist->GetBinError(nx+1,ny+1)));
        }        
    }
    cout << "hist loaded" << endl;
}

void PileupHolder::SetTimeSlice(int n_t_slice) {
    this->n_t_slice = n_t_slice;
}

std::pair<double,double> PileupHolder::GetPDF_t(const int& e1Bin,const int& e2Bin,const int& e3Bin) {
    double p1 = data[n_t_slice-1][e1Bin-1].first;
    double err1 = data[n_t_slice-1][e1Bin-1].second;

    double p2 = data[n_t_slice-1][e2Bin-1].first;
    double err2 = data[n_t_slice-1][e2Bin-1].second;

    double p3 = data[n_t_slice-1][e3Bin-1].first;
    double err3 = data[n_t_slice-1][e3Bin-1].second;

    double val = p1 * p2 * p3 * DT * DT;
    double err = TMath::Sqrt(p1 * p1 * p2 * p2 * err3 * err3
     + p1 * p1 * p3 * p3 * err2 * err2 
     + p2 * p2 * p3 * p3 * err1 * err1) * DT * DT;

    return std::make_pair(val,err);
}

std::pair<double,double> PileupHolder::GetPDF_d(const int& e1Bin,const int& e2Bin) {
    double p1 = data[n_t_slice-1][e1Bin-1].first;
    double err1 = data[n_t_slice-1][e1Bin-1].second;

    double p2 = data[n_t_slice-1][e2Bin-1].first;
    double err2 = data[n_t_slice-1][e2Bin-1].second;

    double val = p1 * p2 * DT;
    double err = TMath::Sqrt(p1 * p1 * err2 * err2 + p2 * p2 * err1 * err1) * DT;
    return std::make_pair(val,err);
}


std::pair<double,double> PileupHolder::GetPDF_t(const double& e_x,const double& e_y,const double& e_z) {
    int e1Bin = FindBin(e_x);
    int e2Bin = FindBin(e_y);
    int e3Bin = FindBin(e_z);
    return GetPDF_t(e1Bin,e2Bin,e3Bin);
}

std::pair<double,double> PileupHolder::GetPDF_d(const double& e_x,const double& e_y) {
    int e1Bin = FindBin(e_x);
    int e2Bin = FindBin(e_y);    
    return GetPDF_d(e1Bin,e2Bin);
}

void ShowUsage(string main_exe) {
    cout << "Usage: \n" << 
        "/Path/to/" << main_exe << " -i input.root -o output.root -s status" << endl;
}

map<string,string> parse_argv(int argc,char *argv[]) {
    map<string,string> config;
    // map<string,string> name_map = {{"-i","input"},{"-o","output"},{"-n","caloNum"},{"-e","energyCut"}};
    map<string,string> name_map = {{"-i","input"},{"-o","output"},{"-s","status"},{"-tstart","tstart"},{"-tend","tend"}};
    int n_arg = 1;
    while(n_arg<argc){
        if(name_map.find(argv[n_arg])!=name_map.end()){
            config[name_map[argv[n_arg]]] = argv[n_arg+1];
            n_arg += 2;
        }
        else {
            ShowUsage(argv[0]);
        }
    }
    if(config.find("tstart")==config.end()) {
        config["tstart"] = "0";
    }
    if(config.find("tend")==config.end()) {
        config["tend"] = std::to_string(nSlices);
    }
    for(auto item : name_map){
        if(config.find(item.second)==config.end()){
            ShowUsage(argv[0]);
            exit(-1);
        }
    }
    return config;
}




TH2F* makePerturbationHist_d(PileupHolder *holder, double maxT, const std::string& name){
    int nslices = holder->GetNSlices();
    TH2F* pert = new TH2F(name.c_str(), name.c_str(), nSlices, 0, maxT, 3 * nPuEnergyBins, 0, 9300);

    auto pertYAxis = pert->GetYaxis();
    auto puE = holder->GetHist()->GetYaxis();
    for(int i_slice = 0; i_slice < nSlices; ++i_slice){
        holder->SetTimeSlice(i_slice+1);
        if(i_slice%100==0) {            
            cout << "Making double perturbation hists - slice = " << i_slice << endl;
        }
        for(int i_bin = 1; i_bin <= 2 * nPuEnergyBins; ++i_bin){
            auto e = pertYAxis->GetBinCenter(i_bin);

            double puSum = 0;
            double puSum_m2 = 0;
            for(int i_e2bin = 1; i_e2bin <= nPuEnergyBins; ++i_e2bin){
                auto e2 = puE->GetBinCenter(i_e2bin);
                auto e1 = e - e2;
                if(e1<0){
                    continue;
                }
                auto i_e1bin = puE->FindBin(e1);
                auto res = holder->GetPDF_d(e1,e2);
                puSum += res.first;
                puSum_m2 += TMath::Power(res.second,2);
            }

            double puDiff = 0;
            double puDiff_m2 = 0;                
            for(int i_e2bin = 1; i_e2bin <= nPuEnergyBins; ++i_e2bin){
                auto e2 = puE->GetBinCenter(i_e2bin);
                auto res = holder->GetPDF_d(e,e2);
                puDiff += res.first;
                puDiff_m2 += TMath::Power(res.second,2);
            }
        pert->SetBinContent(i_slice + 1, i_bin, puSum - 2 * puDiff);
        pert->SetBinError(i_slice + 1, i_bin, TMath::Sqrt(puSum_m2 + 4 * puDiff_m2));
        }
    }
    return pert;
}



TH2F* makePerturbationHist_t(PileupHolder *holder, double maxT, const std::string& name){
    int nslices = holder->GetNSlices();
    TH2F* pert = new TH2F(name.c_str(), name.c_str(), nSlices, 0, maxT, 3 * nPuEnergyBins, 0, 9300);

    auto pertYAxis = pert->GetYaxis();
    auto puE = holder->GetHist()->GetYaxis();
    for(int i_slice = 0; i_slice < nSlices; ++i_slice){
        holder->SetTimeSlice(i_slice+1);
        if(i_slice%100==0) {            
            cout << "Making triple perturbation hists - slice = " << i_slice << endl;
        }
        for(int i_bin = 1; i_bin <= 3*nPuEnergyBins; ++i_bin){
            auto e = pertYAxis->GetBinCenter(i_bin);

            double puSum_s = 0;
            double puSum_s_m2 = 0;
            for(int i_e2bin = 1; i_e2bin <= nPuEnergyBins; ++i_e2bin){
                for(int i_e3bin = 1; i_e3bin <= nPuEnergyBins; ++i_e3bin){
                    auto e2 = puE->GetBinCenter(i_e2bin);
                    auto e3 = puE->GetBinCenter(i_e3bin);
                    auto e1 = e - e2 - e3;
                    if(e1<0) {
                        continue;
                    }

                    auto res = holder->GetPDF_t(e1,e2,e3);                    
                    puSum_s += res.first;
                    puSum_s_m2 += TMath::Power(res.second, 2);
                    
                }
            }

            double puDiff_d = 0;
            double puDiff_d_m2 = 0;
            for(int i_e2bin = 1; i_e2bin <= nPuEnergyBins; ++i_e2bin){
                auto e2 = puE->GetBinCenter(i_e2bin);
                double e1 = e - e2;
                if(e1<0) {
                    continue;
                }
                

                for(int i_e3bin = 1; i_e3bin <= nPuEnergyBins; ++i_e3bin){
                    double e3 = puE->GetBinCenter(i_e3bin);
                    auto res = holder->GetPDF_t(e1,e2,e3);
                    puDiff_d += res.first;
                    puDiff_d_m2 += TMath::Power(res.second, 2);

                }
            }

            double puSum = 0;
            double puSum_m2 = 0;            
            for(int i_e2bin = 1; i_e2bin <= nPuEnergyBins; ++i_e2bin){
                for(int i_e3bin = 1; i_e3bin <= nPuEnergyBins; ++i_e3bin){
                    auto e1 = e;
                    auto e2 = puE->GetBinCenter(i_e2bin);
                    auto e3 = puE->GetBinCenter(i_e3bin);
                    auto res = holder->GetPDF_t(e1,e2,e3);
                    puSum += res.first;
                    puSum_m2 += TMath::Power(res.second, 2);
                }
            }
            pert->SetBinContent(i_slice + 1, i_bin, puSum_s - 3 * puDiff_d + 3 * puSum);
            pert->SetBinError(i_slice + 1, i_bin, TMath::Sqrt(puSum_s_m2 + 9 * puDiff_d_m2 + 9 * puSum_m2));
        }
    }
    return pert;
}

void ScalePileupPerturbation(TH2F * rawHist, TH2 * pert_d, TH2 * pert_t ) {
    double e_start = 3500;
    double e_end = 5000;
    double t_start = 30e3;
    double t_end = 500e3;

    int t_binStart = rawHist->GetXaxis()->FindBin(t_start);
    int t_binEnd = rawHist->GetXaxis()->FindBin(t_end);
    
    
    auto rawE = rawHist->ProjectionY("rawE", t_binStart, t_binEnd);
    auto pertE = pert->ProjectionY("pertE", t_binStart,t_binEnd);
    

    int e_binStart = rawE->GetXaxis()->FindBin(e_start);
    int e_binEnd = rawE->GetXaxis()->FindBin(e_end);
    
    double errRaw = 0;
    double errPert = 0;
    double intRaw = rawE->IntegralAndError(e_binStart,e_binEnd,errRaw);
    double intPert = pertE->IntegralAndError(e_binStart,e_binEnd,errPert);


    double ratio = intRaw/intPert;
    double err_ratio = ratio * TMath::Sqrt(TMath::Power(errRaw/intRaw,2) + TMath::Power(errPert/intPert,2));
    int NX = pert->GetNbinsX();
    int NY = pert->GetNbinsY();
    for(int nx=0;nx<NX;nx++){
        for(int ny=0;ny<NY;ny++){
            double v = pert->GetBinContent(nx+1,ny+1);
            double e = pert->GetBinError(nx+1,ny+1);
            pert->SetBinContent(nx+1,ny+1,v*ratio);
            pert->SetBinError(nx+1,ny+1,TMath::Sqrt(e*e*ratio*ratio+err_ratio*err_ratio*v*v));
        }
    }
}

int main(int argc, char *argv[]){
    auto argv_map = parse_argv(argc,argv);
    
    
    // string status = argv_map["status"];

    std::string hname = "ClusterHitHist/ET_AllCalo1";
    TFile* infile =  TFile::Open(argv_map["input"].c_str());
    cout << argv_map["input"].c_str() << endl;
    

    TH2D* tHist2D = (TH2D*)infile->Get(hname.c_str());
    cout << hname << endl;

    auto nSlices = tHist2D->GetNbinsX();    
    double maxT = tHist2D->GetXaxis()->GetBinLowEdge(nSlices + 1);

    TH2F* perturbedDist = (TH2F*) tHist2D->Clone();
    TH2F* lastRho = perturbedDist;
    perturbedDist->SetName("rhoc");
    perturbedDist->SetTitle("rhoc");
 
    string slice_label = Form("Slice%s_%s",argv_map["tstart"],argv_map["tend"]);
    string outf_name = Form("%s_%s.root",argv_map["output"],slice_label)
    TFile* outf = new TFile(outf_name.c_str(), "recreate");
    outf->cd();
    perturbedDist->Write();

    PileupHolder * holder = new PileupHolder();
    holder->SetETHist(lastRho);

    
    // triple correction
    // up to three iterations
    for(int i = 0; i < 1; i++){
        auto empPert_d = makePerturbationHist_d(holder, maxT, Form("deltaRho%i_d", i+1));
        auto empPert_t = makePerturbationHist_t(holder, maxT, Form("deltaRho%i_t", i+1));

        
        ScalePileupPerturbation(perturbedDist, empPert_d, empTrueDist_t);

        TH2F* empTrueDist_d = new TH2F(Form("rhoc%i_t", i+1), Form("rhoc%i_t", i+1), nSlices, 0, maxT, 3 * nPuEnergyBins, 0, 9300);

        for(int j = 1; j <= empTrueDist_d->GetNbinsX(); j++){
            for(int k = 1; k <= 2 * nPuEnergyBins; k++){
                empTrueDist_d->SetBinContent(j, k, perturbedDist->GetBinContent(j, k) - empPert_d->GetBinContent(j, k));

                double errEmp_m2 = TMath::Power(perturbedDist->GetBinError(j, k), 2) + TMath::Power(empPert_d->GetBinError(j, k), 2);
                empTrueDist_d->SetBinError(j, k, TMath::Sqrt(errEmp_m2));
            }
        }
        outf->cd();

        empPert_d->Write();
        empPert_t->Write();
        empTrueDist_d->Write();

        lastRho = empTrueDist_d;
    }

    outf->Close();
    return 0; 
}