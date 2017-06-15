#!/usr/bin/env python

#extend left flanking
def extend_left_mis(seq, start, left, mlen, motif, offset1, offset2):
	i = left + offset1
	while i >= 0:
		j = ((i - start + offset2) % mlen + mlen) % mlen
		if seq[i] != motif[j]:
			return i + 1
		i -= 1
	return i

def extend_left_sub(seq, start, left, mlen, motif):
	return extend_left_mis(seq, start, left, mlen, motif, -1, 0)

def extend_left_ins(seq, start, left, mlen, motif):
	return extend_left_mis(seq, start, left, mlen, motif, -1, 1)

def extend_left_del(seq, start, left, mlen, motif):
	return extend_left_mis(seq, start, left, mlen, motif, 0, -1)

#extend right flank
def extend_right_mis(seq, seqlen, start, right, mlen, motif, offset1, offset2):
	i = right + offset1
	while i < seqlen:
		j = (i - start + offset2) % mlen
		if seq[i] != motif[j]:
			return i - 1
		i += 1
	return i

def extend_right_sub(seq, seqlen, start, right, mlen, motif):
	return extend_right_mis(seq, seqlen, start, right, mlen, motif, 1, 0)

def extend_right_ins(seq, seqlen, start, right, mlen, motif):
	return extend_right_mis(seq, seqlen, start, right, mlen, motif, 1, -1)

def extend_right_del(seq, seqlen, start, right, mlen, motif):
	return extend_right_mis(seq, seqlen, start, right, mlen, motif, 0, 1)


sequence = 'AAGAAGAAGACGAGACAGAAGAAG'

total_len = len(sequence)

i = 0
seed_start = 0
seed_length = 0

while i < total_len:
	j = 1
	while j < 7:
		seed_start = i
		seed_length = j

		while seed_start+seed_length < total_len and sequence[i] == sequence[i+j]:
			i += 1
			seed_length += 1

		repeat = seed_length/j

		if repeat >= 3:
			motif = sequence[seed_start:seed_start+j]
			
			#the start and end position of imperfect ssr
			start = seed_start
			end = seed_start + seed_length - 1

			errors = 0 # consecutive errors
			
			matches = seed_length
			substitution = 0
			deletion = 0
			insertion = 0
			
			sub_index = 0
			ins_index = 0
			del_index = 0

			#extend left flank sequence
			extend_start = seed_start
			extend_index = seed_start - 1
			while 1:
				if extend_index < 0:
					break

				errors += 1
				sub_index = extend_left_sub(sequence, extend_start, extend_index, j, motif)
				ins_index = extend_left_ins(sequence, extend_start, extend_index, j, motif)
				del_index = extend_left_del(sequence, extend_start, extend_index, j, motif)

				sub_extend = extend_index - sub_index
				ins_extend = extend_index - ins_index
				del_extend = extend_index - del_index

				min_index = min((sub_index, ins_index, del_index))
				max_extend = max((sub_extend, ins_extend, del_extend))

				if max_extend > 0:
					errors = 0
					matches += max_extend
					start = min_index
					if max_extend == sub_extend:
						substitution += 1
					elif max_extend == ins_extend:
						insertion += 1
					else:
						deletion += 1

				if errors > 3:
					break

				if min_index == sub_index:
					extend_index = sub_index - 1
				elif min_index == ins_index:
					extend_index = ins_index - 1
					extend_start -= 1
				else:
					extend_index = del_index - 1
					extend_start += 1

			#extend right flank
			errors = 0
			extend_start = seed_start
			extend_index = seed_start + seed_length
			while 1:
				if extend_index >= total_len:
					break

				errors += 1
				sub_index = extend_right_sub(sequence, total_len, extend_start, extend_index, j, motif)
				ins_index = extend_right_ins(sequence, total_len, extend_start, extend_index, j, motif)
				del_index = extend_right_del(sequence, total_len, extend_start, extend_index, j, motif)

				sub_extend = sub_index - extend_index
				ins_extend = ins_index - extend_index
				del_extend = del_index - extend_index

				max_index = max((sub_index, ins_index, del_index))
				max_extend = max((sub_extend, ins_extend, del_extend))

				if max_extend > 0:
					errors = 0
					matches += max_extend
					end = max_index
					if max_extend == sub_extend:
						substitution += 1
					elif max_extend == ins_extend:
						insertion += 1
					else:
						deletion += 1

				if errors > 3:
					break

				if max_index == sub_index:
					extend_index = sub_index + 1
				elif max_index == ins_index:
					extend_index = ins_index + 1
					extend_start += 1
				else:
					extend_index = del_index + 1
					extend_start -= 1

			gap = insertion + deletion
			if gap:
				score = matches - substitution - 5 - (gap - 1)*2
			else:
				score = matches - substitution

			if score >= 5:
				print motif, j, start+1, end+1, end-start+1, matches, substitution, insertion, deletion, score
				i = end + 1
			else:
				i = seed_start + 1

			j = 0
		else:
			i = seed_start
			
		j += 1

	i += 1


