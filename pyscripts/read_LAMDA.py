from astropy.table import Table
from astropy import units as u


def read_LAMDA_file(fle,store_collisions=None,store_headers=None):
    '''
    Reads in any LAMDA-structured file, stores each table as an astropy table and returns
    the collection of tables in a dictionary. Two optional flags will activate 
    the storage of collisional rates and/or headers for reprinting back to file.

    Strucure of LAMDA file described here: https://home.strw.leidenuniv.nl/~moldata/molformat.html

    Collisions may be with multiple partners, so collision tables are stored in another dictionary
    of tables, with each key corresponding to LAMDA-valid identifications:
        1=H2, 2=para-H2, 3=ortho-H2, 4=electrons, 5=H, 6=He, 7=H+

    output: dict{'levels','transitions','collisions','headers'}

    '''
    output_dict = {}
    header_dict = {}
    try:
        testopen = open(fle,'r')
        testopen.close()
    except IOError:
        print("File not found!")
    else:
        with open(fle,'r') as f:
            level_header = [f.readline() for _ in range(7)]
            header_dict['levels'] = level_header
            # Store number of levels for table
            num_levels = int(level_header[5])

            # Determine columns present in level table for this molecule
            level_column_labels = ['level','energy','g_ul']
            num_level_columns =  len(level_header[6].split('+'))
            # Extra columns may include quantum numbers to describe molecular energy level
            if num_level_columns > 3:
                level_column_labels += [label.strip() for label in level_header[6].split('+')[3:]]

            # Assume quantum number columns tabulated in integer format
            level_column_dtypes = [int,float,float]+[int]*(num_level_columns-3)
            level_column_dtypes_quantum_str = [int,float,float]+[str]*(num_level_columns-3)

            quantum_string_flag = 0
            apply_dtypes = lambda x,y: x(y) 
            table_values = [f.readline().split() for _ in range(num_levels)]

            # Attempt to type each data entry with dtype defined above
            # If failure present, presume quantum number issue and leave as str
            try:
                for n_row in range(len(table_values)):
                    if (num_level_columns==len(table_values[n_row])):
                        table_values[n_row] = [apply_dtypes(x,y) for x,y in zip(level_column_dtypes,table_values[n_row])]
                    else:
                        raise Exception("Row {0} in the level table has an incorrect number of entries!".format(n_row))
            except ValueError:
                print("Default formatting of level table failed: leaving quantum numbers as strings\n")
                quantum_string_flag = 1
                for n_row in range(len(table_values)):
                    table_values[n_row] = [apply_dtypes(x,y) for x,y in \
                                            zip(level_column_dtypes_quantum_str,table_values[n_row])]
            else:
                level_table = Table(rows=table_values,names=level_column_labels,dtype=level_column_dtypes)

            if quantum_string_flag:
                level_table = Table(rows=table_values,names=level_column_labels,dtype=level_column_dtypes_quantum_str)

            level_table['energy'] /= u.cm

            output_dict['levels'] = level_table


            transition_header = [f.readline() for _ in range(3)]
            header_dict['transitions'] = transition_header
            num_transitions = int(transition_header[1])

            # Determine columns present in transition table - default four
            transition_column_labels = ['transition','upper','lower','A_ul']
            num_transition_columns = len(transition_header[2].split('+'))
            # Extra columns often include frequency in GHz, E_u in K
            if num_transition_columns > 4:
                transition_column_labels += [label.strip() for label in transition_header[2].split('+')[4:]]

            ## Construct transitions table

            # Assume extra columns, if any, can be typed to float
            transition_column_dtypes = [int,int,int,float]+[float]*(num_transition_columns-4)
            transition_column_dtypes_str = [int,int,int,float]+[str]*(num_transition_columns-4)
            transition_string_flag = 0

            table_values = [f.readline().split() for _ in range(num_transitions)]

            # Attempt to type each data entry with dtype defined above
            # If failure present, leave as strings
            try:
                for n_row in range(len(table_values)):
                    if (num_transition_columns==len(table_values[n_row])):
                        table_values[n_row] = [apply_dtypes(x,y) for x,y in zip(transition_column_dtypes,table_values[n_row])]
                    else:
                        raise Exception("Row {0} in the transition table has an incorrect number of entries!".format(n_row))
            except ValueError:
                print("Default formatting of transition table failed: leaving extra columns as strings\n")
                transition_string_flag = 1
                for n_row in range(len(table_values)):
                    table_values[n_row] = [apply_dtypes(x,y) for x,y in \
                                            zip(transition_column_dtypes_str,table_values[n_row])]
            else:
                transition_table = Table(rows=table_values,names=transition_column_labels,dtype=transition_column_dtypes)

            if transition_string_flag:
                transition_table = Table(rows=table_values,names=transition_column_labels,
                                          dtype=transition_column_dtypes_str)

            transition_table['A_ul'] /= u.s

            output_dict['transitions'] = transition_table


            # Determine if we need to read in collision tables
            # This can be very time-consuming, so it's worth avoiding if possible.
            if store_collisions:
                # Create dictionary to provide the variable number of collision tables in a single package.
                collisions_dict = {}
                collision_block_header = [f.readline() for _ in range(2)]
                header_dict['collision_block'] = collision_block_header

                num_collision_partners = int(collision_block_header[1])

                # For each partner present in file, extract collision table
                for n_partner in range(num_collision_partners):
                    # Header provides partner id, number of collisions and number of columns(temperatures)
                    collision_partner_header = [f.readline() for _ in range(9)]
                    partner_id = int(collision_partner_header[1].split()[0])
                    header_dict['collision_partner{0}'.format(partner_id)] = collision_partner_header
                    num_collisions = int(collision_partner_header[3])
                    num_collision_temperatures = int(collision_partner_header[5])

                    # Number of columns in table is three plus number of temperatures given
                    num_collision_columns = num_collision_temperatures+3

                    # To label columns, need to apply temperatures given
                    collision_column_labels = ['transition','upper','lower']+['k_{0}'.format(x) for x in \
                                                collision_partner_header[7].split()]

                    collision_column_dtypes = [int,int,int]+[float]*num_collision_temperatures

                    table_values = [f.readline().split() for _ in range(num_collisions)]


                    for n_row in range(len(table_values)):
                        if (num_collision_columns==len(table_values[n_row])):
                            table_values[n_row] = [apply_dtypes(x,y) for x,y in \
                                                    zip(collision_column_dtypes,table_values[n_row])]
                        else:
                            print(table_values[n_row])
                            raise Exception("Row {0} in the transition table has an incorrect number of entries!".format(n_row))
            
                    # Store into collision table
                    collision_table = Table(rows=table_values,names=collision_column_labels,dtype=collision_column_dtypes)

                    # Assign units to rates in table
                    for column in collision_table.colnames:
                        if 'k' in column:
                            collision_table[column] *= (u.cm**3/u.s)

                    collisions_dict[partner_id] = collision_table

                output_dict['collisions'] = collisions_dict

        if store_headers:
            output_dict['headers'] = header_dict
            
        return output_dict
