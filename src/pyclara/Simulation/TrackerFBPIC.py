import pyclara

def elegant_lte_splitter(filename, split_element, output_file1, output_file2):

    # Step 1: load the full file
    elements = pyclara.Converters.elegant_lte_loader(filename)

    # Step 2: find LINE key and the ordered element list
    line_key  = next(k for k, v in elements.items() if v['TYPE'] == 'LINE')
    full_line = elements[line_key]['LINE']

    # Step 3: find split index — exact match
    if split_element not in full_line:
    
    # check if a case-insensitive match exists
        close_match = next(
            (e for e in full_line if e.upper() == split_element.upper()), 
            None
        )
    
        if close_match:
            raise ValueError(
                f"'{split_element}' not found. Did you mean '{close_match}'?"
            )
        else:
            raise ValueError(
            f"'{split_element}' not found in beamline."
            )

    split_idx = full_line.index(split_element)


    line1 = full_line[:split_idx]    # before split element
    line2 = ['START'] + full_line[split_idx:]    # from split element onwards

    # Step 4: build two element dicts
    set1 = set(e.upper() for e in line1)
    set2 = set(e.upper() for e in line2)

    elements1 = {}
    elements2 = {}

    elements2 = {'START': elements['START']} 

    for name, elem in elements.items():
        if elem['TYPE'] == 'LINE':
            continue
        if name.upper() in set1:
            elements1[name] = elem
        elif name.upper() in set2:
            elements2[name] = elem


    # Step 5: add LINE entries to each dict
    elements1[line_key + '_1'] = {'NAME': line_key + '_1', 'TYPE': 'LINE', 'LINE': line1}
    elements2[line_key + '_2'] = {'NAME': line_key + '_2', 'TYPE': 'LINE', 'LINE': line2}

    # Step 6: write both files
    pyclara.Converters.elegant_lte_writer(elements1, output_file1)
    pyclara.Converters.elegant_lte_writer(elements2, output_file2)



