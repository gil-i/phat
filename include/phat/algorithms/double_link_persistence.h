/*  Copyright 2024 Purdue University
    Contributed by: Gilberto Gonzalez

    This file is NOT part of PHAT.

    THIS FILE IMPLEMENTS A NEW STATE-OF-THE-ART ALGORITHM to compute 
    the double-linked persistence of a system of a filtration, 
    it is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    If not, see <http://www.gnu.org/licenses/>. */

#pragma once

#include <phat/helpers/misc.h>
#include <phat/boundary_matrix.h>
#include <iostream>

template <typename S>
std::ostream& operator<<(std::ostream& os,
                    const std::vector<S>& vector)
{
    // Printing all the elements
    // using <<
    for (auto element : vector) {
        os << element << " ";
    }
    return os;
}


namespace phat {
    class double_linked_reduction { public:
        template< typename Representation >
        void operator() ( boundary_matrix< Representation >& boundary_matrix ) {

            const index nr_columns = boundary_matrix.get_num_cols();
            std::vector< index > lowest_one_lookup( nr_columns, -1 );
            phat::boundary_matrix <Representation> V_matrix; // The matrix V, in the decomposition of the reduced matrix has the form R=DV
            V_matrix.set_num_cols( nr_columns );

            for( index cur_col = 0; cur_col < nr_columns; cur_col++ ) {
                index lowest_one = boundary_matrix.get_max_index( cur_col );
                column col_temp;
                while( lowest_one != -1 && lowest_one_lookup[ lowest_one ] != -1 ) {
                    boundary_matrix.add_to( lowest_one_lookup[ lowest_one ], cur_col );
                    V_matrix.add_to(lowest_one_lookup[ lowest_one], cur_col);
                    col_temp.push_back(lowest_one_lookup[ lowest_one ]);
                    //std::cout << "cur_col="<<cur_col<<", lowest_one_lookup="<<lowest_one_lookup[ lowest_one]<< ", lowest_one="<<lowest_one <<std::endl;
                    lowest_one = boundary_matrix.get_max_index( cur_col );
                }
                // adding col_temp to cur_col
                std::sort( col_temp.begin(), col_temp.end() );
                V_matrix.set_col(0 ,col_temp);
                V_matrix.add_to(0  , cur_col);
                V_matrix.clear(0); // resetting first column to 0 (it should always be zero)
                std::vector<index> col_print(nr_columns);
                std::cout<< col_print <<std::endl;
                V_matrix.get_col(cur_col, col_print);
                std::cout << "Cur col ="<< cur_col << std::endl;
                std::cout << col_temp << std::endl;
                std::cout<< col_print <<std::endl;
                if( lowest_one != -1 ) {
                    lowest_one_lookup[ lowest_one ] = cur_col;
                }
            }
            
            // clearing paired columns
            for( index idx = 0; idx < boundary_matrix.get_num_cols(); idx++ ) { 
                if( !boundary_matrix.is_empty( idx ) ) {
                    V_matrix.clear(idx);
                }
            }
            boundary_matrix = V_matrix;
            standard_reduction red2;
            red2(boundary_matrix);
            std::cout << "boundary_matrix" << std::endl;
        }
    };
}

