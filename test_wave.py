import wave
import numpy
from wave import TWOPI
import matplotlib.pyplot as plt


LMBD = 0.4
KVEC = numpy.array([TWOPI/LMBD, 0., 0.])
POINT = numpy.array([-1.2, 0., 0.])


def test_gradIncident():
	nvec = numpy.array([-1., 0., 0.])
	res = wave.gradIncident(nvec, KVEC, POINT)
	print('test_gradIncident: {}'.format(res))

def test_computeScatteredWaveElement():
	p0 = numpy.array([1.,  0.1, 0.])
	p1 = numpy.array([0.9,-0.1, 0.])
	res = wave.computeScatteredWaveElement(KVEC, p0, p1, POINT)
	print('test_computeScatteredWaveElement: {}'.format(res))

def test_computeScatteredWave():
	nt = 100
	nt1 = nt + 1
	dt = TWOPI / float(nt)
	xc = numpy.array([numpy.cos(i*dt) for i in range(nt1)])
	yc = numpy.array([numpy.sin(i*dt) for i in range(nt1)])
	res = wave.computeScatteredWave(KVEC, xc, yc, POINT)
	print('test_computeScatteredWave: {}'.format(res))

def test_computeScatteredWave_front():
	nt = 100
	nt1 = nt + 1
	dt = TWOPI / float(nt)
	xc = numpy.array([numpy.cos(i*dt) for i in range(nt1)])
	yc = numpy.array([numpy.sin(i*dt) for i in range(nt1)])
	no = 40
	xo = [-1.05-0.02*i for i in range(no)]
	inci = [wave.incident(KVEC, numpy.array([xo[i], 0., 0.])) for i in range(no)]
	inci_real = numpy.array([s.real for s in inci])
	inci_imag = numpy.array([s.imag for s in inci])
	scat = [wave.computeScatteredWave(KVEC, xc, yc, numpy.array([xo[i], 0., 0.])) for i in range(no)]
	scat_real = numpy.array([s.real for s in scat])
	scat_imag = numpy.array([s.imag for s in scat])
	plt.plot(xo, scat_real, 'b-.')
	plt.plot(xo, scat_imag, 'r-.')
	plt.plot(xo, inci_real, 'b--')
	plt.plot(xo, inci_imag, 'r--')
	plt.plot(xo, inci_real + scat_real, 'b-')
	plt.plot(xo, inci_imag + scat_imag, 'r-')

	plt.legend(['sct re', 'sct im', 'inc re', 'inc im', 'tot re', 'tot im'])
	plt.xlabel('x observer')
	plt.title('front')

def test_computeScatteredWave_back():
	nt = 100
	nt1 = nt + 1
	dt = TWOPI / float(nt)
	xc = numpy.array([numpy.cos(i*dt) for i in range(nt1)])
	yc = numpy.array([numpy.sin(i*dt) for i in range(nt1)])
	no = 40
	xo = [1.05+0.02*i for i in range(no)]
	inci = [wave.incident(KVEC, numpy.array([xo[i], 0., 0.])) for i in range(no)]
	inci_real = numpy.array([s.real for s in inci])
	inci_imag = numpy.array([s.imag for s in inci])
	scat = [wave.computeScatteredWave(KVEC, xc, yc, numpy.array([xo[i], 0., 0.])) for i in range(no)]
	scat_real = numpy.array([s.real for s in scat])
	scat_imag = numpy.array([s.imag for s in scat])
	plt.plot(xo, scat_real, 'b-.')
	plt.plot(xo, scat_imag, 'r-.')
	plt.plot(xo, inci_real, 'b--')
	plt.plot(xo, inci_imag, 'r--')
	plt.plot(xo, inci_real + scat_real, 'b-')
	plt.plot(xo, inci_imag + scat_imag, 'r-')

	plt.legend(['sct re', 'sct im', 'inc re', 'inc im', 'tot re', 'tot im'])
	plt.xlabel('x observer')
	plt.title('back')

def test_box():
	nt = 100
	dx, dy = 2./float(nt), 2./float(nt)
	xc = [1. for i in range(nt)] + [1. - i*dx for i in range(nt)] + [-1. for i in range(nt)] + [-1. + i*dx for i in range(nt + 1)]
	yc = [-1. + i*dy for i in range(nt)] + [1. for i in range(nt)] + [1. - i*dy for i in range(nt)] + [-1. for i in range(nt + 1)]
	xc - numpy.array(xc)
	yc = numpy.array(yc)
	xoL = numpy.linspace(-2., -1.01, 40)
	xoR = numpy.linspace(1.01, 2., 40)
	noL, noR = len(xoL), len(xoR)

	inciL = [wave.incident(KVEC, numpy.array([xoL[i], 0.1, 0.])) for i in range(noL)]
	inciL_real = numpy.array([s.real for s in inciL])
	inciL_imag = numpy.array([s.imag for s in inciL])

	inciR = [wave.incident(KVEC, numpy.array([xoR[i], 0.1, 0.])) for i in range(noR)]
	inciR_real = numpy.array([s.real for s in inciR])
	inciR_imag = numpy.array([s.imag for s in inciR])

	scatL = [wave.computeScatteredWave(KVEC, xc, yc, numpy.array([xoL[i], 0.1, 0.])) for i in range(noL)]
	scatL_real = numpy.array([s.real for s in scatL])
	scatL_imag = numpy.array([s.imag for s in scatL])

	scatR = [wave.computeScatteredWave(KVEC, xc, yc, numpy.array([xoR[i], 0.1, 0.])) for i in range(noR)]
	scatR_real = numpy.array([s.real for s in scatR])
	scatR_imag = numpy.array([s.imag for s in scatR])


	plt.plot(xoL, scatL_real, 'b-.')
	plt.plot(xoL, scatL_imag, 'r-.')
	plt.plot(xoL, inciL_real, 'b--')
	plt.plot(xoL, inciL_imag, 'r--')
	plt.legend(['sct re', 'sct im', 'inc re', 'inc im', 'tot re', 'tot im'])
	plt.plot(xoR, scatR_real, 'b-.')
	plt.plot(xoR, scatR_imag, 'r-.')
	plt.plot(xoR, inciR_real, 'b--')
	plt.plot(xoR, inciR_imag, 'r--')

	plt.xlabel('x observer')
	plt.title('box')



if __name__ == '__main__':
	test_gradIncident()
	test_computeScatteredWaveElement()
	test_computeScatteredWave()
	test_computeScatteredWave_front()
	test_computeScatteredWave_back()
	test_box()
	plt.show()
