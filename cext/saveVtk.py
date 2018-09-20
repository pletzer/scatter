
def saveData(filename, x, y, data, fieldname='data'):
	ny1, nx1 = data.shape
	with open(filename, 'w') as f:
		f.write('# vtk DataFile Version 3.0\n')
		f.write('vtk output\nASCII\n')
		f.write('DATASET RECTILINEAR_GRID\n')
		f.write('DIMENSIONS {} {} 1\n'.format(nx1, ny1))
		f.write('X_COORDINATES {} double\n'.format(nx1))
		for i in range(nx1):
			f.write('{} '.format(x[i]))
		f.write('\n')
		f.write('Y_COORDINATES {} double\n'.format(ny1))
		for j in range(ny1):
			f.write('{} '.format(y[j]))
		f.write('\n')
		f.write('Z_COORDINATES 1 double\n0\n')
		f.write('POINT_DATA {}\n'.format(ny1*nx1))
		f.write('FIELD FieldData 1\n')
		f.write('{} 1 {} float\n'.format(fieldname, ny1*nx1))
		for j in range(ny1):
			for i in range(nx1):
				f.write('{} '.format(data[j, i]))
			f.write('\n')
		f.close()