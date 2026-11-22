#pragma once
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
namespace panmanUtils {
// Progress belongs on stderr so it never becomes FASTA/MAF data. Rate-limit
// intermediate updates, but always show phase boundaries and completion.
class ConstructionProgress {
    using Clock = std::chrono::steady_clock;
    Clock::time_point started = Clock::now(), previous = started;
public:
    void update(const std::string& phase, size_t done, size_t total,
                const std::string& detail = {}, bool force = false) {
        auto now = Clock::now();
        if(!force && now-previous < std::chrono::seconds(5)) return;
        previous = now;
        std::ostringstream message;
        message << "[MAF] " << phase << ": " << done;
        if(total) message << '/' << total << " (" << std::fixed << std::setprecision(1)
                          << 100.0 * double(done)/double(total) << "%)";
        if(!detail.empty()) message << "; " << detail;
        message << "; elapsed " << std::chrono::duration_cast<std::chrono::seconds>(now-started).count() << "s";
        std::cerr << message.str() << std::endl;
    }
};
}
