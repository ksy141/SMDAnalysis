NC3 = ['N',   'C12',  'H12A', 'H12B',
       'C13', 'H13A', 'H13B', 'H13C',
       'C14', 'H14A', 'H14B', 'H14C',
       'C15', 'H15A', 'H15B', 'H15C',
       'C11', 'H11A', 'H11B']

PO4 = ['P', 'O13',  'O14', 'O12', 'O11']

GL  = ['C1', 'HA', 'HB',
       'C2', 'HS', 
       'O21', 'C21', 'O22',
       'C3',  'HX',  'HY',
       'O31', 'C31', 'O32']

C21  = ['C2%d' %n for n in range(2, 7)]
C21 += ['H%dR' %n for n in range(2, 7)]
C21 += ['H%dS' %n for n in range(2, 7)]

C22  = ['C2%d' %n for n in range(7, 11)]
C22 += ['H%dR' %n for n in range(7, 11)]
C22 += ['H%dS' %n for n in range(7, 11)]
C22 += ['H91', 'H101']

C23  = ['C2%d' %n for n in range(11, 15)]
C23 += ['H%dR' %n for n in range(11, 15)]
C23 += ['H%dS' %n for n in range(11, 15)]

C24  = ['C2%d' %n for n in range(15, 19)]
C24 += ['H%dR' %n for n in range(15, 19)]
C24 += ['H%dS' %n for n in range(15, 19)]
C24 += ['H18T']

C31  = ['C3%d' %n for n in range(2, 7)]
C31 += ['H%dX' %n for n in range(2, 7)]
C31 += ['H%dY' %n for n in range(2, 7)]

C32  = ['C3%d' %n for n in range(7, 11)]
C32 += ['H%dX' %n for n in range(7, 11)]
C32 += ['H%dY' %n for n in range(7, 11)]

C33  = ['C3%d' %n for n in range(11, 15)]
C33 += ['H%dX' %n for n in range(11, 15)]
C33 += ['H%dY' %n for n in range(11, 15)]

C34  = ['C3%d' %n for n in range(15, 17)]
C34 += ['H%dX' %n for n in range(15, 17)]
C34 += ['H%dY' %n for n in range(15, 17)]
C34 += ['H16Z']

group = {'NC3': NC3,
         'PO4': PO4,
         'GL':  GL,
         'C21': C21,
         'C22': C22,
         'C23': C23,
         'C24': C24,
         'C31': C31,
         'C32': C32,
         'C33': C33,
         'C34': C34}
