#include <vector>

// Forward declaration of CellState (assumed already defined)
class CellState {
    public:
        int Ncell;      // Number of cells
        int Nface;      // Number of faces (typically equal to Ncell for periodic boundaries)
        double dx;      // Cell width
        double gamma;   // Ratio of specific heats
    
        // Conserved variables, one value per cell.
        std::vector<double> rho;
        std::vector<double> momentum;
        std::vector<double> energy;
    
        // Cell face arrays: each face stores two values:
        // index 0: state from the cell to the left (right face of cell i)
        // index 1: state from the cell to the right (left face of cell i+1)
        std::vector< std::vector<double> > rho_face;
        std::vector< std::vector<double> > momentum_face;
        std::vector< std::vector<double> > energy_face;
    
        // Flux arrays at faces
        std::vector<double> flux_rho;
        std::vector<double> flux_momentum;
        std::vector<double> flux_energy;
    
        // Example constructor (details omitted; see previous answer)
        CellState(int Ncell_, double dx_,
                  const std::vector<double>& rho_,
                  const std::vector<double>& momentum_,
                  const std::vector<double>& energy_,
                  double gamma_ = 5.0/3.0)
           : Ncell(Ncell_), dx(dx_), gamma(gamma_),
             rho(rho_), momentum(momentum_), energy(energy_)
        {
             _allocate_cell_faces();
        }
    
    private:
        void _allocate_cell_faces() {
             // For periodic boundaries in 1D, we set Nface equal to Ncell.
             Nface = Ncell;
             rho_face.resize(Nface, std::vector<double>(2, 0.0));
             momentum_face.resize(Nface, std::vector<double>(2, 0.0));
             energy_face.resize(Nface, std::vector<double>(2, 0.0));
    
             flux_rho.resize(Nface, 0.0);
             flux_momentum.resize(Nface, 0.0);
             flux_energy.resize(Nface, 0.0);
        }
    };