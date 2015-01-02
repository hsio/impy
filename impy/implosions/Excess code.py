
    # it turns out there is NO need to parse the zone information section
    # however, I already written the code for it, so it is just commented out
    # in case it is wanted in the future

                # # if we run into the zoning information, record it
                # elif self.lilacKeys["KEY_ZONING_HEADER"] in line:
                #     self.lilac_full_parse_zone_info(dataReader)
                #     self.runProgress += 0.1



    # def lilac_full_parse_zone_info(self,dataReader):
    #     """Parse the zoning information block."""
    #     # if this function is called, assume it is at beginning of zoning info block
    #
    #     #first, fetch the next line,  and time should be the first
    #     while True:
    #         try:
    #             line = next(dataReader)
    #             line = line[0]
    #         except:
    #             continue
    #
    #         if line.split() == []:
    #             continue
    #
    #         # if first word is index, we arrived at the headings
    #         elif self.lilacKeys["KEY_ZONE_DATA_INDEX"] == line.split()[0]:
    #             # now we need to find the order of the headings
    #             headings = self.split_str_into_words_min_two_spaces(line)
    #             print("\nheading for zone info")
    #             print(headings,"\n")
    #
    #         # hey, are you a  information row?
    #         elif line.split()[0].isdigit():
    #             line = line.split()
    #             line = [self.float_(item) for item in line]
    #             # attempt to insert data into respective array
    #             try:
    #                 temp_index = headings.index(self.lilacKeys["KEY_ZONE_INFO_IONS"])
    #                 num_ni = line[temp_index]
    #                 temp_index = headings.index(self.lilacKeys["KEY_ZONE_INFO_VOL"])
    #                 vol = line[temp_index]
    #                 self.niRaw.append(num_ni / vol)
    #             except:
    #                 print( "something wrong in parsing zoning information")
    #                 continue
    #
    #         elif self.lilacKeys["KEY_ZONING_FOOTER"] in line:
    #             zoning_info_processed = True
    #             break