def print_step_header(title):
	print('------ [' + title + ']')

def fn(x):
	"""adds double quotes before and after the passed string, to sanitize filenames"""
	return ('"' + x + '"')
