#include <iostream>
#include <vector>
#include <fstream> // Include for file handling
#include <iomanip> // Include for fixed precision formatting

using namespace std;

const int N = 30; // The size of the matrix (excluding boundary)

const double length = 10;

double h = length / N;

const double dt = 1;
const double diff = 0.1;

const double source = 100;
double a = diff * dt / (h * h);

vector<vector<double>> u(N + 2, vector<double>(N + 2));
vector<vector<double>> u0(N, vector<double>(N));

vector<vector<double>> v(N + 2, vector<double>(N + 2));
vector<vector<double>> v0(N + 2, vector<double>(N + 2));

vector<vector<double>> dens(N + 2, vector<double>(N + 2));
vector<vector<double>> dens0(N + 2, vector<double>(N + 2));

void SWAP(vector<vector<double>> &dens, vector<vector<double>> &dens0)
{
    vector<vector<double>> temp(N + 2, vector<double>(N + 2));
    temp = dens;
    dens0 = temp;
}

void addSource(vector<vector<double>> &dens)
{
    for (int i = 20; i <= 25; i++)
    {
        for (int j = 20; j <= 25; j++)
        {
            dens0[i][j] = source * dt;
        }
    }

    for (int i = 3; i <= 8; i++)
    {
        for (int j = 3; j <= 8; j++)
        {
            dens0[i][j] = source * dt;
        }
    }
}

void diffuse(vector<vector<double>> &dens, vector<vector<double>> &dens0)
{

    for (int k = 0; k < 300; k++)
    {
        for (int i = 1; i <= N; i++)
        {
            for (int j = 1; j <= N; j++)
            {
                dens[i][j] = (dens0[i][j] + a * (dens[i - 1][j] + dens[i + 1][j] + dens[i][j - 1] + dens[i][j + 1])) / (1 + 4 * a);
            }
        }
    }
}

// Function to save dens0 to a file as a matrix in CSV format for Excel
// Function to save dens0 to a file as a matrix in CSV format for Excel
void saveToFile(const vector<vector<double>> &dens, const string &filename)
{
    ofstream outFile(filename);
    if (outFile.is_open())
    {
        // Set fixed point notation and set precision for writing to the file
        outFile << fixed << setprecision(6); // Set the precision to 6 decimal places

        // Write the data row by row, each row being a line in the CSV
        for (int i = 0; i <= N + 1; i++) // Include boundary cells (0 to N+1)
        {
            for (int j = 0; j <= N + 1; j++) // Include boundary cells (0 to N+1)
            {
                outFile << dens[i][j]; // Write the value

                if (j < N + 1)      // Avoid adding a comma at the end of the row
                    outFile << ","; // Separate values with a comma
            }
            outFile << "\n"; // New line for each row
        }

        outFile.close();
        cout << "Data saved to " << filename << endl;
    }
    else
    {
        cerr << "Unable to open file: " << filename << endl;
    }
}

int main()
{

    for (int t = 0; t < 100; t = t + dt)
    {
        addSource(dens);

        saveToFile(dens0, "dens0_t" + to_string(t) + ".csv");
        // Print dens0 in the console for each timestep (optional)
        for (int i = 0; i <= N + 1; i++) // Include boundary cells (0 to N+1)
        {
            for (int j = 0; j <= N + 1; j++) // Include boundary cells (0 to N+1)
            {
                cout << dens0[i][j] << " ";
            }
            cout << "\n";
        }

        diffuse(dens, dens0);
        SWAP(dens, dens0);
        cout << "\n";
        cout << "\n";

        // Save the current dens0 to a file after each time step
    }
    cout << a;

    return 0;
}
