/*! Merge several files into one filtering duplicated entries.
*/

#include "TauMLTools/Core/interface/program_main.h"
#include "TauMLTools/Analysis/interface/TauTuple.h"
#include "TauMLTools/Analysis/interface/L2TrainingTauTuple.h"
#include "TauMLTools/Analysis/interface/L2EventTuple.h"
#include "TauMLTools/Core/interface/RootFilesMerger.h"

struct Arguments {
    run::Argument<std::string> output{"output", "output root file"};
    run::Argument<std::vector<std::string>> input_dirs{"input-dir", "input directory"};
    run::Argument<std::string> file_name_pattern{"file-name-pattern", "regex expression to match file names",
                                                 "^.*\\.root$"};
    run::Argument<std::string> exclude_list{"exclude-list", "comma separated list of files to exclude", ""};
    run::Argument<std::string> exclude_dir_list{"exclude-dir-list",
                                                "comma separated list of directories to exclude", ""};
    run::Argument<unsigned> n_threads{"n-threads", "number of threads", 1};
};

namespace analysis {

class MergeL2Tuples : public RootFilesMerger {
public:
    using TrainingTau = tau_train_tuple::Tau;
    using TrainingTauTuple = tau_train_tuple::TauTrainTuple;
    using EntryId = tau_train_tuple::L2TrainTupleEntryId;
    using EntryIdSet = std::set<EntryId>;

    MergeL2TauTuples(const Arguments& args) :
        RootFilesMerger(args.output(), args.input_dirs(), args.file_name_pattern(), args.exclude_list(),
                        args.exclude_dir_list(), args.n_threads(), ROOT::kZLIB, 9),
        output_tauTuple("taus", output_file.get(), false), //output_summaryTuple("summary", output_file.get(), false),
        n_total_duplicates(0)
    {
    }

    void Run()
    {
        Process(false, false);

        output_tauTuple.Write();

        std::cout << "All file has been merged. Number of files = " << input_files.size()
                  << ". Number of output entries = " << output_tauTuple.GetEntries()
                  << ". Total number of duplicated entires = " << n_total_duplicates << "." << std::endl;
    }

private:
    virtual void ProcessFile(const std::string& /*file_name*/, const std::shared_ptr<TFile>& file) override
    {
        TrainingTauTuple input_tauTuple("taus", file.get(), true);
        size_t n_duplicates = 0;
        for(const TrainingTau& tau : input_tauTuple) {
            const EntryId entry_id(tau);
            if(processed_entries.count(entry_id)) {
                ++n_duplicates;
                continue;
            }
            processed_entries.insert(entry_id);
            output_tauTuple() = tau;
            output_tauTuple.Fill();
        }
        n_total_duplicates += n_duplicates;
        /*
        SummaryTuple input_summaryTuple("summary", file.get(), true);
        for(const ProdSummary& summary : input_summaryTuple) {
            output_summaryTuple() = summary;
            output_summaryTuple.Fill();
        }
        */
        std::cout << "\tn_entries = " << input_tauTuple.GetEntries() << ", n_duplicates = " << n_duplicates << ".\n";
    }

private:
    TrainingTauTuple output_tauTuple; 
    EntryIdSet processed_entries;
    size_t n_total_duplicates;
};

} // namespace analysis

PROGRAM_MAIN(analysis::MergeL2Tuples, Arguments)
