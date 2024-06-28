#include <vector>
#include <iostream>
#include <bitset>
#include <iomanip>
#include <random>
#include <cassert>

using namespace std;

uint64_t to_binary(double x) {
    return * (uint64_t*) (&x);
}

double exponential_inv(double cdf) {
    // cdf = 1 - exp(-x)
    return -std::log(1 - cdf);
}

std::vector<bool> encode(double center, double lambda, double v) {
    // For now assume center = 0 and lambda = 1

    double left = 0;
    double right = 1;
    uint64_t diff = std::numeric_limits<uint64_t>::max();

    std::vector<bool> output;

    for (int i = 0; i < 64; ++i) {
        const double right_val = exponential_inv(right);
        const double left_val = exponential_inv(left);
        uint64_t new_diff = to_binary(right_val) - to_binary(left_val);
        const double ratio = new_diff / double(diff);
        if (new_diff >= diff || new_diff < 8) {

            // should return a simple binary encoding of the left over
            assert(to_binary(v) >= to_binary(left_val));

            const bitset<64> simple_diff = (to_binary(v) - to_binary(left_val));
            uint64_t extra_bits = 0;
            uint64_t cutoff = 1;
            while (cutoff < new_diff) {
                extra_bits += 1;
                cutoff *= 2;
            }

            std::cout << "simple_diff: " << simple_diff << std::endl;

            for (uint64_t bit_idx = 0; bit_idx < extra_bits; ++bit_idx) {
                std::cout << "simple: " << bit_idx << " " << simple_diff[bit_idx] << std::endl;
                output.push_back(simple_diff[bit_idx]); // TODO check endian
            }
            std::cout << "diff didn't change! " << new_diff << " " << diff << std::endl;

            return output;
        } else {
            diff = new_diff;
        }





        const double mid = (left + right) / 2;

        const double guess = exponential_inv(mid);
        if (v > guess) {
            left = mid;
            output.push_back(false);
        } else {
            right = mid;
            output.push_back(true);
        }


        std::cout << i << " " << guess << " " << v << " " << right_val - left_val << " " << to_binary(right_val) - to_binary(left_val) << " " << ratio << std::endl;
    }

}

int main() {
    std::cout << std::setprecision(14);


    std::random_device rd; 
    std::mt19937 gen(rd()); 

    //default_random_engine generator;
    exponential_distribution<double> distribution(0.1);
    //normal_distribution<double> distribution(0.0, 1);
    static constexpr size_t nbits = 64;

    std::vector<int> count(nbits, 0.0);

    double a = distribution(gen);
    const auto output = encode(0, 1, a);
    std::cout << "got size: " << output.size() << std::endl;
}