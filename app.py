import os
import datetime
from flask import Flask, request, session, render_template, redirect, url_for, send_from_directory, flash
from werkzeug import secure_filename
import phylomito
import shortuuid

app = Flask(__name__,static_url_path='/data/')

app.config['PREFERRED_URL_SCHEME'] = 'http'
app.config['SECRET_KEY'] = os.environ['SECRET_KEY']
app.config['UPLOAD_FOLDER'] = '/data/uploads/'
app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024
ALLOWED_EXTENSIONS = set(['gbk', 'gb'])
IGNORED_FILES = set(['.gitignore'])
default_args = {'log': 'phyml.log', 'extension': ['.gbk', '.gb'], 'bootstrap': 100, 'protein': False, 'gene_tree': False, 'dloop': False, 'code_table': 2}

def allowed_file(filename):
    print filename.rsplit('.', 1)[1].lower(), ALLOWED_EXTENSIONS
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def gen_file_name(filename):
    """
    If file was exist already, rename it and return a new name
    """

    i = 1
    while os.path.exists(os.path.join(app.config['UPLOAD_FOLDER'], filename)):
        name, extension = os.path.splitext(filename)
        filename = '%s_%s%s' % (name, str(i), extension)
        i += 1

    return filename

@app.template_filter()
def datetimefilter(value, format='%Y/%m/%d %H:%M'):
    """convert a datetime to a different format."""
    return value.strftime(format)

app.jinja_env.filters['datetimefilter'] = datetimefilter

@app.route("/", methods=['GET', 'POST'])
def template_root():	
    if request.method == 'POST':
        if 'file' not in request.files:
            return redirect(request.url)
        files = []
        for i in request.files.items(True):
            file = i[1]
                # if user does not select file, browser also
            # submit a empty part without filename
            if file.filename == '':
                return redirect(request.url)
            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                files.append(filename)
                file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            else:
                print 'file not allowed'
            if len(files) > 1:
                session['job_id'] = shortuuid.uuid()
                outpath = app.config['UPLOAD_FOLDER'] + session['job_id']+ '/'
                os.mkdir(outpath)
                args = default_args.copy()
                args['inpath'] = app.config['UPLOAD_FOLDER']
                args['outpath'] = outpath
                phylomito.main(args)
    return render_template('welcome.html') 

#                uploaded_file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
#                size = os.path.getsize(uploaded_file_path)
#        files = [f for f in os.listdir(app.config['UPLOAD_FOLDER']) if os.path.isfile(os.path.join(app.config['UPLOAD_FOLDER'],f)) and f not in IGNORED_FILES ]
#            size = os.path.getsize(os.path.join(app.config['UPLOAD_FOLDER'], f))
#

@app.route("/delete/<string:filename>", methods=['DELETE'])
def delete(filename):
    file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
    if os.path.exists(file_path):
        try:
            os.remove(file_path)
        except:
            pass
            

if __name__ == '__main__':
    app.run('127.0.0.1', 8000, debug=True)


