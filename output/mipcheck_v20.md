# MIP signal vs real beam dots in single_eval_v20_dz4

Steep tracks (dip > 60 deg), 80 events; per track the median darkness (max within +-2 px, gray levels) and the persistence (fraction of slices darker than 6); cells: 16 / 50 / 84 %. Real dots are selected by persistence >= 0.6; 'passing' applies the same cut to the signal; 'followed' measures the signal like the real dots (followed from its position in the middle slice).

| sample | n | darkness | persistence | passing | darkness (passing) | followed: darkness, persistence |
|---|---|---|---|---|---|---|
| real persistent dots | 128082 | 22.6 / 40.6 / 55.3 | 0.906 / 1 / 1 | - | - | - |
| followed random positions | 16000 | 10.3 / 21.8 / 40.9 | 0.615 / 0.882 / 1 | 85% | - | - |
| fixed random positions | 16000 | 5.49 / 11.3 / 19.5 | 0.48 / 0.706 / 0.902 | 68% | - | - |
| signal dE/dx 0-0.0005 MeV/um | 55 | 10.8 / 14.2 / 24.8 | 0.739 / 0.846 / 0.964 | 100% | 10.8 / 14.2 / 24.8 | 13 / 24.1 / 44.9, 0.709 / 0.914 / 1 |
| signal dE/dx 0.0005-0.0015 MeV/um | 25 | 15.1 / 20.5 / 32.5 | 0.853 / 0.93 / 1 | 100% | 15.1 / 20.5 / 32.5 | 18.4 / 26.1 / 59.3, 0.801 / 1 / 1 |
| signal dE/dx 0.0015-1 MeV/um | 111 | 31.8 / 66.4 / 123 | 1 / 1 / 1 | 100% | 31.8 / 66.4 / 123 | 22.4 / 34.2 / 60.1, 0.853 / 1 / 1 |
