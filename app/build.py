src = open('app.src.html', encoding='utf-8').read()
model = open('model.js', encoding='utf-8').read()
assert '/*MODEL*/' in src
out = src.replace('/*MODEL*/', model)
open('../index.html', 'w', encoding='utf-8').write(out)
print(len(out), 'bytes')
