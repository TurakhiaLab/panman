#include <iostream>
#include <fstream>
#include <string>
#include <cctype>

int main(int argc, char ** argv) {
    std::ifstream inputFile(argv[1]);
    std::ofstream outputFile(argv[2]);

    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open the input file." << std::endl;
        return 1;
    }

    if (!outputFile.is_open()) {
        std::cerr << "Error: Could not open the output file." << std::endl;
        return 1;
    }

    char ch;
    while (inputFile.get(ch)) { // Reads each character from the input file
        if (std::islower(ch)) {
            ch = std::toupper(ch); // Capitalize the character if it's lowercase
        }
        outputFile.put(ch); // Write the character to the output file
    }

    inputFile.close();
    outputFile.close();

    std::cout << "File processing completed. Check output.txt for results." << std::endl;

    return 0;
}
