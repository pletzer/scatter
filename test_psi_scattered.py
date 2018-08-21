import psi_scattered
import numpy
from psi_scattered import TWOPI


LMBD = 0.4
KVEC = numpy.array([TWOPI/LMBD, 0., 0.])
POINT = numpy.array([-1.2, 0., 0.])


def test_gradPsiIncident():
	nvec = numpy.array([-1., 0., 0.])
	res = psi_scattered.gradPsiIncident(nvec, KVEC, POINT)
	print('test_gradPsiIncident: {}'.format(res))

def test_computePsiScatteredElement():
	p0 = numpy.array([1., 0., 0.])
	p1 = numpy.array([0.9, 0.5, 0.])
	res = psi_scattered.computePsiScatteredElement(KVEC, p0, p1, POINT)
	print('test_computePsiScatteredElement: {}'.format(res))

def test_computePsiScattered():
	nt = 10
	nt1 = nt + 1
	dt = TWOPI / float(nt)
	xc = numpy.array([numpy.cos(i*dt) for i in range(nt1)])
	yc = numpy.array([numpy.sin(i*dt) for i in range(nt1)])
	res = psi_scattered.computePsiScattered(KVEC, xc, yc, POINT)
	print('test_computePsiScattered: {}'.format(res))

if __name__ == '__main__':
	test_gradPsiIncident()
	test_computePsiScatteredElement()
	test_computePsiScattered()