class Renumerator(object):
    """
    Generic renumerator class. Provides methods for changing the ID numbering
    in all ModernaStructure-based objects.
    """
    def __init__(self, struct):
        """
        :Arguments:
        * struct - the structure to be renumbered (descendant of ModernaStructure)
        """
        #struct = validate_structure(struct)
        self.struct = struct
        self.original_ids = [i.id[1] for i in self.struct]
        self.ids = self.original_ids

    def apply_alignment(self, alignment):
        """
        Inserts gaps basing on the input alignment.
        Please note that the SECOND sequence is used.
        """
        seq = str(alignment.aligned_sequences[1])
        seq += "\xff"  # a little hack to avoid "open" rightmost gaps
        gaps = []
        gap_length = 0
        shift = 0	# remember that each insertion shifts the numbering!
        res_num = None
        res_iter = self.ids.__iter__()
        for s in seq:
            if s=="-":
                if gap_length==0:
                    start_pos = res_num
                gap_length += 1
            elif s!="_":
                if s!="\xff":
                    res_num = next(res_iter)
                if gap_length > 0:
                    gaps.append((start_pos + shift, gap_length))
                    shift += gap_length
                    gap_length = 0
        for gap in gaps:
            self.insert_gap(gap[0], gap[1])

    def insert_gap(self, position, length):
        """
        Inserts a gap of a given length at a specified position.
        """
        if position>=len(self.ids):
            return
        elif position==0:
            gap_id = 1
            new_ids = []
        else:
            gap_id = self.ids[position - 1] + 1
            new_ids = self.ids[:position]
        shift = length
        for i in self.ids[position:]:
            if i>gap_id:
                shift -= i - gap_id
                gap_id = i
                if shift<0: shift = 0
            new_ids.append(i + shift)
            gap_id += 1
        self.ids = new_ids

    def renumber(self, count_from):
        """
        Provides simple renumeration of the entire structure,
        starting from count_from with an increment of 1.
        """
        id_gen = self.id_generator(count_from)
        self.ids = [next(id_gen) for i in self.ids]

    def id_generator(self, first_id):
        """
        Dummy ID generator, needs to be replaced with a working code.
        """
        pass

    def apply(self):
        """
        Applies all changes to the input structure.
        """
        resi_list = list(self.struct)
        for r in self.struct:
            self.struct.remove_residue(r.identifier)
        for r, i in zip(resi_list, self.ids):
            self.struct.add_residue(r, str(i))


class NumberRenumerator(Renumerator):
    """
    Renumerator based on number-only IDs.
    """
    def id_generator(self, first_id):
        i = first_id
        while True:
            yield i
            i += 1

