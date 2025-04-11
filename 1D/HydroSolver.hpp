#include <vector>
#include <cassert>

class CellState;

class HydroSolver {
    public:
        // Pointer to the cell state and simulation time variables
        CellState* cellstate;
        double t0;
        double t;
    
        // Constructor: store a pointer to the CellState and initialize time
        HydroSolver(CellState* cellstate_, double t0_ = 0.0)
           : cellstate(cellstate_), t0(t0_), t(t0_) { }
    
        // Run the simulation until time t1 with timestep dt
        void run(double t1, double dt) {
            assert(t1 > t);
            while (t < t1) {
                double dt_local = std::min(dt, t1 - t);
                _step(dt_local);
            }
        }
    
    private:
        // Single time step: reconstruction, flux computation, application, and time update
        void _step(double dt) {
            _reconstruction();
            _compute_fluxes();
            _apply_fluxes(dt);
            t += dt;
        }
    
        // Reconstruction: for each face, assign the left state from cell i and the right state from cell (i+1)
        void _reconstruction() {
            for (int i = 0; i < cellstate->Nface; i++) {
                cellstate->rho_face[i][0]      = cellstate->rho[i];
                cellstate->momentum_face[i][0] = cellstate->momentum[i];
                cellstate->energy_face[i][0]   = cellstate->energy[i];
    
                // For periodic boundaries, use modulo arithmetic to get cell i+1.
                int ip1 = (i + 1) % cellstate->Ncell;
                cellstate->rho_face[i][1]      = cellstate->rho[ip1];
                cellstate->momentum_face[i][1] = cellstate->momentum[ip1];
                cellstate->energy_face[i][1]   = cellstate->energy[ip1];
            }
        }

        // Apply the computed fluxes to update the cell conserved variables.
        void _apply_fluxes(double dt) {
            double dx = cellstate->dx;
            for (int i = 0; i < cellstate->Nface; i++) {
                int ip1 = (i + 1) % cellstate->Ncell;
                cellstate->rho[i]      -= (dt / dx) * cellstate->flux_rho[i];
                cellstate->rho[ip1]    += (dt / dx) * cellstate->flux_rho[i];
                cellstate->momentum[i] -= (dt / dx) * cellstate->flux_momentum[i];
                cellstate->momentum[ip1] += (dt / dx) * cellstate->flux_momentum[i];
                cellstate->energy[i]   -= (dt / dx) * cellstate->flux_energy[i];
                cellstate->energy[ip1] += (dt / dx) * cellstate->flux_energy[i];
            }
        }
    
        // Compute fluxes at each face using the Rusanov (local Lax-Friedrichs) solver.
        void _compute_fluxes() {
            for (int i = 0; i < cellstate->Nface; i++) {
                // Extract left and right states at the face
                double rho_L      = cellstate->rho_face[i][0];
                double rho_R      = cellstate->rho_face[i][1];
                double momentum_L = cellstate->momentum_face[i][0];
                double momentum_R = cellstate->momentum_face[i][1];
                double energy_L   = cellstate->energy_face[i][0];
                double energy_R   = cellstate->energy_face[i][1];
    
                double flux_rho, flux_momentum, flux_energy;
                _rusanov_flux(rho_L, rho_R, momentum_L, momentum_R, energy_L, energy_R,
                              flux_rho, flux_momentum, flux_energy);
    
                cellstate->flux_rho[i]      = flux_rho;
                cellstate->flux_momentum[i] = flux_momentum;
                cellstate->flux_energy[i]   = flux_energy;
            }
        }
    
        // Helper: Compute the Euler flux for a given state.
        // The fluxes are returned via reference arguments.
        void _Euler_flux(double rho, double momentum, double energy,
                         double &flux_rho, double &flux_momentum, double &flux_energy) {
            double vel = momentum / rho;
            double kinetic_energy = 0.5 * rho * vel * vel;
            double pressure = (energy - kinetic_energy) * (cellstate->gamma - 1.0);
            flux_rho = rho * vel;
            flux_momentum = momentum * vel + pressure;
            flux_energy = (energy + pressure) * vel;
        }
    
        // Helper: Compute the Rusanov flux at a face given left and right states.
        void _rusanov_flux(double rho_L, double rho_R, double momentum_L, double momentum_R,
                           double energy_L, double energy_R,
                           double &flux_rho, double &flux_momentum, double &flux_energy) {
            double flux_rho_L, flux_momentum_L, flux_energy_L;
            double flux_rho_R, flux_momentum_R, flux_energy_R;
            _Euler_flux(rho_L, momentum_L, energy_L, flux_rho_L, flux_momentum_L, flux_energy_L);
            _Euler_flux(rho_R, momentum_R, energy_R, flux_rho_R, flux_momentum_R, flux_energy_R);
    
            // Compute sound speeds for the face states.
            double vel_L = momentum_L / rho_L;
            double kinetic_energy_L = 0.5 * rho_L * vel_L * vel_L;
            double pressure_L = (energy_L - kinetic_energy_L) * (cellstate->gamma - 1.0);
            double cs_L = std::sqrt(cellstate->gamma * pressure_L / rho_L);
    
            double vel_R = momentum_R / rho_R;
            double kinetic_energy_R = 0.5 * rho_R * vel_R * vel_R;
            double pressure_R = (energy_R - kinetic_energy_R) * (cellstate->gamma - 1.0);
            double cs_R = std::sqrt(cellstate->gamma * pressure_R / rho_R);
    
            // Maximum wave speed (alpha) at the interface.
            double alpha = std::max(std::abs(vel_L) + cs_L, std::abs(vel_R) + cs_R);
    
            flux_rho      = 0.5 * (flux_rho_L + flux_rho_R) - 0.5 * alpha * (rho_R - rho_L);
            flux_momentum = 0.5 * (flux_momentum_L + flux_momentum_R) - 0.5 * alpha * (momentum_R - momentum_L);
            flux_energy   = 0.5 * (flux_energy_L + flux_energy_R) - 0.5 * alpha * (energy_R - energy_L);
        }
    };