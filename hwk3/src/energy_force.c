#include "energy_force.h"
#include "params.h"
#include "atoms.h"
#include "timer.h"
#include <mpi.h>

//************************************************************************
// compute_long_range_correction() function
//   - Calculates long range correction due to finite interaction cutoff.
//   - Arguments:
//       - len_jo: struct containing leonard jones interaction parameters.
//       - m_pars: struct containing misc. parameters.
//       - energy_long: long-range correction to energy.
//       - force_long: long-range correction to force.
//************************************************************************
void compute_long_range_correction(const lj_params * len_jo, const misc_params * m_pars,
                                   float * energy_long, float * force_long )
{

   float ulongpre = m_pars->float_N * 8.0 * len_jo->eps * 
                    m_pars->pi * m_pars->density;
   *energy_long = ulongpre * ( len_jo->sig12 / ( 9.0 * len_jo->rcut9 ) - 
                               len_jo->sig6 / ( 6.0 * len_jo->rcut3 ) );

   float vlongpre = 96.0 * len_jo->eps * m_pars->pi * m_pars->density;
   *force_long = -1.0 * vlongpre * ( len_jo->sig12 / ( 9.0 * len_jo->rcut9 ) - 
                                     len_jo->sig6 / ( 6.0 * len_jo->rcut3 ) ); 

}

//************************************************************************
// compute_energy_and_force() function
//   - Calculates energy and force acting on each atom.
//   - Arguments:
//       - myatoms: struct containing all atomic information.
//       - len_jo: struct containing lennard jones interaction parameters.
//       - m_pars: struct containing misc. parameters.
//************************************************************************
void compute_energy_and_force( Atoms * myatoms, const lj_params * len_jo, 
                               const misc_params * m_pars )
{

   timeit(1,0);

   int atomi, atomj;
   float pot_energy[myatoms->N];
   float virial[myatoms->N];
   #pragma simd
   for (atomi=0; atomi < myatoms->N; atomi++)
   {
      myatoms->fx[atomi] = 0.0f;
      myatoms->fy[atomi] = 0.0f;
      myatoms->fz[atomi] = 0.0f;
   }
   myatoms->pot_energy = 0.0;
   myatoms->virial = 0.0;

   MPI_Allgatherv(MPI_IN_PLACE,myatoms->N,MPI_FLOAT,myatoms->xx,myatoms->recvcounts,myatoms->displs,MPI_FLOAT,MPI_COMM_WORLD);
   MPI_Allgatherv(MPI_IN_PLACE,myatoms->N,MPI_FLOAT,myatoms->yy,myatoms->recvcounts,myatoms->displs,MPI_FLOAT,MPI_COMM_WORLD);
   MPI_Allgatherv(MPI_IN_PLACE,myatoms->N,MPI_FLOAT,myatoms->zz,myatoms->recvcounts,myatoms->displs,MPI_FLOAT,MPI_COMM_WORLD);

   for (atomi=myatoms->start_index; atomi < myatoms->end_index; atomi++)
   {

      #pragma simd
      for (atomj=0 ; atomj < myatoms->N; atomj++)
      {

         pot_energy[atomj] = 0.0f;
         virial[atomj] = 0.0f;

         if ( atomi != atomj )
         {

            float xxi = myatoms->xx[atomi] - myatoms->xx[atomj];
            xxi = minimum_image( xxi, m_pars->side, m_pars->sideh );
            float yyi = myatoms->yy[atomi] - myatoms->yy[atomj];
            yyi = minimum_image( yyi, m_pars->side, m_pars->sideh );
            float zzi = myatoms->zz[atomi] - myatoms->zz[atomj];
            zzi = minimum_image( zzi, m_pars->side, m_pars->sideh );

            float dis2 = xxi*xxi + yyi*yyi + zzi*zzi;
            if ( dis2 <= len_jo->rcut2 )
            {

               float dis2i = 1.0 / dis2;
               float dis6i = dis2i * dis2i * dis2i;
               float dis12i = dis6i * dis6i;
               pot_energy[atomj] = len_jo->sig12 * dis12i - 
                                   len_jo->sig6 * dis6i;
               float fterm = dis2i * ( 2.0 * len_jo->sig12 * dis12i -
                                          len_jo->sig6 * dis6i );
               virial[atomj] = fterm * dis2;
            
               myatoms->fx[atomj] -= fterm * xxi;
               myatoms->fy[atomj] -= fterm * yyi;
               myatoms->fz[atomj] -= fterm * zzi;
         
            }
    
         }

      } 

      for (atomj=0 ; atomj < myatoms->N; atomj++)
      {
         myatoms->pot_energy += pot_energy[atomj];
         myatoms->virial -= virial[atomj];
      }

   }

   #pragma simd
   for (atomi=0; atomi < myatoms->N; atomi++)
   {
      myatoms->fx[atomi] *= 24.0 * len_jo->eps;
      myatoms->fy[atomi] *= 24.0 * len_jo->eps;
      myatoms->fz[atomi] *= 24.0 * len_jo->eps;
   }

   myatoms->pot_energy *= 2.0 * len_jo->eps; // factor of 1/2 included here
   myatoms->virial *= 12.0 * len_jo->eps; // factor of 1/2 included here

   timeit(1,1);

}

//**********************************************************************
// minimum_image() function
//   - Finds the nearest images of atom i and j, and returns distance.
//   - Arguments:
//       - dist: 1d distance between atoms i and j in central sim. cell.
//       - box_length: length of simulation cell.
//       - half_box_length: half of length of simulation cell.
//**********************************************************************
float minimum_image( const float dist, const float box_length, 
                     const float half_box_length )
{

   float min_dist = dist;
   if (dist > half_box_length ) min_dist = dist - box_length; 
   if (dist < -half_box_length ) min_dist = dist + box_length;
   return min_dist; 

}