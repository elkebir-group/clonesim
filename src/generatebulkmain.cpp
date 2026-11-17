/*
 * generatebulk: draw bulk sequencing data from the output of `simulate`.
 *
 * Reads node.tsv and proportions.tsv of a simulated instance and writes
 *   snvs.tsv  snv, sample, segment, variant, total
 *   cn.tsv    segment, sample, x, y, fraction
 *
 * In sample p, the mean copy number of segment s is c = sum_v (x_v + y_v) u_v,
 * where u_v is the proportion of node v in p. The total reads of an SNV on s are
 * Poisson with mean depth * c / 2 (so a diploid segment has the given depth),
 * and its variant reads are binomial with probability sum_v (xbar_v + ybar_v) u_v / c.
 * The copy-number calls are the fraction of cells in each state (x, y) of each
 * segment. With -cn_noise > 0, each fraction gets Gaussian noise of that standard
 * deviation, and the fractions are clipped at 0 and renormalized.
 */
#include "utils.h"
#include "arg_parser.h"
#include <boost/random.hpp>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::vector<std::map<std::string, std::string>> read_tsv(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("cannot open " + filename);
    std::string line;
    std::vector<std::string> header;
    std::vector<std::map<std::string, std::string>> rows;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::vector<std::string> fields;
        std::stringstream ss(line);
        std::string f;
        while (std::getline(ss, f, '\t')) fields.push_back(f);
        if (header.empty()) {
            header = fields;
            continue;
        }
        if (fields.size() != header.size()) throw std::runtime_error("wrong number of fields in " + filename + ": " + line);
        std::map<std::string, std::string> row;
        for (size_t i = 0; i < header.size(); ++i) row[header[i]] = fields[i];
        rows.push_back(row);
    }
    return rows;
}

}  // namespace

int main(int argc, char** argv) {
    std::string in_dir, out_dir;
    double depth = 100;
    double cn_noise = 0;
    int seed = 0;

    lemon2::ArgParser ap(argc, argv);
    ap.refOption("in_dir", "Directory with the output of simulate (node.tsv, proportions.tsv)", in_dir, true)
            .refOption("out_dir", "Output directory for snvs.tsv and cn.tsv (default: in_dir)", out_dir, false)
            .refOption("depth", "Mean total reads of an SNV on a diploid segment (default: 100)", depth, false)
            .refOption("cn_noise", "Standard deviation of Gaussian noise on the copy-number state fractions (default: 0)", cn_noise, false)
            .refOption("s", "Random number generator seed (default: 0)", seed, false);
    ap.parse();
    if (out_dir.empty()) out_dir = in_dir;

    try {
        if (depth <= 0) throw std::runtime_error("-depth must be positive");
        if (cn_noise < 0) throw std::runtime_error("-cn_noise must be nonnegative");
        g_rng = boost::random::mt19937(seed);

        // proportions: node, then one column per sample
        std::ifstream pin(in_dir + "/proportions.tsv");
        if (!pin) throw std::runtime_error("cannot open " + in_dir + "/proportions.tsv");
        std::string line;
        std::getline(pin, line);
        std::map<int, std::vector<double>> prop;
        int nr_samples = -1;
        while (std::getline(pin, line)) {
            if (line.empty()) continue;
            std::stringstream ss(line);
            int node;
            ss >> node;
            std::vector<double> u;
            double x;
            while (ss >> x) u.push_back(x);
            if (nr_samples < 0) nr_samples = u.size();
            if ((int) u.size() != nr_samples || nr_samples == 0)
                throw std::runtime_error("proportions.tsv: every node needs one proportion per sample");
            prop[node] = u;
        }

        // states of every node and segment, and mutated copies of every node and SNV
        std::map<std::pair<int, int>, std::pair<int, int>> state;  // (node, segment) -> (x, y)
        std::map<int, int> snv_segment;
        std::map<std::pair<int, int>, int> mutated;                 // (node, snv) -> xbar + ybar
        for (const auto& r : read_tsv(in_dir + "/node.tsv")) {
            int v = std::stoi(r.at("node")), s = std::stoi(r.at("segment")), m = std::stoi(r.at("m"));
            if (!prop.count(v)) throw std::runtime_error("node " + std::to_string(v) + " has no proportions");
            state[{v, s}] = {std::stoi(r.at("x")), std::stoi(r.at("y"))};
            if (m >= 0) {
                snv_segment[m] = s;
                mutated[{v, m}] = std::stoi(r.at("xbar")) + std::stoi(r.at("ybar"));
            }
        }

        std::map<int, std::vector<double>> mean_cn;  // segment -> per sample
        std::map<int, std::vector<std::map<std::pair<int, int>, double>>> fraction;
        for (const auto& [vs, xy] : state) {
            auto [v, s] = vs;
            auto& c = mean_cn[s];
            auto& f = fraction[s];
            c.resize(nr_samples, 0.0);
            f.resize(nr_samples);
            for (int p = 0; p < nr_samples; ++p) {
                c[p] += (xy.first + xy.second) * prop[v][p];
                f[p][xy] += prop[v][p];
            }
        }

        std::ofstream snvs(out_dir + "/snvs.tsv");
        if (!snvs) throw std::runtime_error("cannot write " + out_dir + "/snvs.tsv");
        snvs << "snv\tsample\tsegment\tvariant\ttotal\n";
        for (const auto& [m, s] : snv_segment) {
            for (int p = 0; p < nr_samples; ++p) {
                double y = 0;
                for (const auto& [v, u] : prop) {
                    auto it = mutated.find({v, m});
                    if (it != mutated.end()) y += it->second * u[p];
                }
                double c = mean_cn[s][p];
                int total = 0, variant = 0;
                if (c > 0) {
                    boost::random::poisson_distribution<int> reads(depth * c / 2);
                    total = reads(g_rng);
                    boost::random::binomial_distribution<int> var(total, std::min(1.0, std::max(0.0, y / c)));
                    variant = var(g_rng);
                }
                snvs << m << '\t' << p << '\t' << s << '\t' << variant << '\t' << total << '\n';
            }
        }

        std::ofstream cn(out_dir + "/cn.tsv");
        if (!cn) throw std::runtime_error("cannot write " + out_dir + "/cn.tsv");
        cn << "segment\tsample\tx\ty\tfraction\n";
        boost::random::normal_distribution<double> noise(0.0, cn_noise > 0 ? cn_noise : 1.0);
        for (auto& [s, per_sample] : fraction) {
            for (int p = 0; p < nr_samples; ++p) {
                auto& f = per_sample[p];
                if (cn_noise > 0) {
                    double sum = 0;
                    for (auto& [xy, x] : f) {
                        x = std::max(0.0, x + noise(g_rng));
                        sum += x;
                    }
                    for (auto& [xy, x] : f) x = sum > 0 ? x / sum : 1.0 / f.size();
                }
                for (const auto& [xy, x] : f)
                    cn << s << '\t' << p << '\t' << xy.first << '\t' << xy.second << '\t' << x << '\n';
            }
        }
        std::cerr << "Wrote " << snv_segment.size() << " SNVs and " << fraction.size() << " segments in "
                  << nr_samples << " samples to " << out_dir << std::endl;
    }
    catch (std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
