import os
import pandas as pd

from flask import Flask, render_template, request, jsonify, json
from werkzeug.utils import secure_filename, redirect

# import dash
# from dash.dependencies import Input, Output
# import dash_table
# import dash_core_components as dcc
# import dash_html_components as html

import genedataset as gd
import compoundmodel as cm
# import dash_plottable as dash

dir = os.path.abspath(os.getcwd())
dir = (dir+'\\file\\')
dir = dir.replace('\\', '/')
server = Flask(__name__)
server.config['UPLOAD_FOLDER'] = dir
server.add_url_rule('/file/<path:filename>', endpoint='file', view_func=server.send_static_file)

filename = ""


@server.route('/')
def index():
    return render_template('/index.html')


@server.route('/dataset', methods=['GET', 'POST'])
def dataset():
    if request.method == 'POST':

        filename = ''
        if request.form.get('hiddenfilename2') is None:
            print('upload dataset')
            file = request.files['upload']
            filename = secure_filename(file.filename)
        else:
            print('sample dataset')
            filename = request.form.get('hiddenfilename2')
        fileexistance = gd.check_exist(filename)

        html = ""

        if fileexistance == True:
            df_headerdataset = gd.read_headerdataset(filename)
            print(df_headerdataset)
            innerhtml = """<form class ="form" action="http://localhost:5000/result" method="POST" enctype="multipart/form-data">"""
            innerhtml += """<table style="text-align: center;"><tr><th>Sample Name</th><th>Sample Type</th></tr>"""
            i = 0
            for header in list(df_headerdataset):
                if i != 0:
                    innerhtml += """<tr><td>""" + header + """</td><td>"""
                    innerhtml += """<input id="toggle-on""" + str(
                        i) + """" class="toggle toggle-left" name="toggle""" + str(
                        i) + """" value="Normal" type="radio" checked>"""
                    innerhtml += """<label for="toggle-on""" + str(i) + """" class="btnradio">Normal</label>"""
                    innerhtml += """<input id="toggle-off""" + str(
                        i) + """" class="toggle toggle-right" name="toggle""" + str(
                        i) + """" value="Perturbation" type="radio">"""
                    innerhtml += """<label for="toggle-off""" + str(
                        i) + """" class="btnradio">Perturbation</label></td></tr>"""
                i += 1
            innerhtml += """</table><br>\
            <div class="modal fade" id="exampleModal" tabindex="-1" role="dialog" aria-labelledby="exampleModalLabel" aria-hidden="true">\
                  <div class="modal-dialog" role="document">\
                    <div class="modal-content">\
                      <div class="modal-header">\
                        <h5 class="modal-title" id="exampleModalLabel">File has already existed</h5>\
                        <button type="button" class="close" data-dismiss="modal" aria-label="Close">\
                          <span aria-hidden="true">&times;</span>\
                        </button>\
                      </div>\
                      <div class="modal-body model_left">\
                          <p>This file has already existed in the system, would you like to update the setting?<p/>\
                            <div class="form-check">\
                              <input class="form-check-input" type="radio" name="settinguser" id="exampleRadios1" value="no" checked>\
                              <label class="form-check-label" for="exampleRadios1">\
                                Continue with the previous setting\
                              </label>\
                            </div>\
                            <div class="form-check">\
                              <input class="form-check-input" type="radio" name="settinguser" id="exampleRadios2" value="yes">\
                              <label class="form-check-label" for="exampleRadios2">\
                                Update new setting\
                              </label>\
                            </div>\
                      </div>\
                      <div class="modal-footer">\
                        <button type="button" class="btn3" data-dismiss="modal">Cancel</button>\
                        <button type="submit" class="btn" name="submit" id="analyze" value="analyze">Submit</button>\
                      </div>\
                    </div>\
                  </div>\
                </div>\
            """
            innerhtml += """<button type="button" class="btn" data-toggle="modal" data-target="#exampleModal" id="btntrigger"><i class="fas fa-search"></i> Analyze</button>"""
            innerhtml += """<input type="hidden" id="hiddenfilename" name="hiddenfilename" value=\"""" + filename + """\"></form>"""
            html = """
            <script type="text/javascript">
            """
            html += """document.getElementById('datasetconfig').innerHTML = '""" + innerhtml + """';"""
            html += """
            </script>
            """
        else:
            file.save(os.path.join(server.config['UPLOAD_FOLDER'], filename))
            df_headerdataset = gd.read_headerdataset(filename)
            print(df_headerdataset)
            innerhtml = """<form class ="form" action="http://localhost:5000/result" method="POST" enctype="multipart/form-data">"""
            innerhtml += """<table style="text-align: center;"><tr><th>Sample Name</th><th>Sample Type</th></tr>"""
            i = 0
            for header in list(df_headerdataset):
                if i != 0:
                    innerhtml += """<tr><td>""" + header + """</td><td>"""
                    innerhtml += """<input id="toggle-on""" + str(
                        i) + """" class="toggle toggle-left" name="toggle""" + str(
                        i) + """" value="Normal" type="radio" checked>"""
                    innerhtml += """<label for="toggle-on""" + str(i) + """" class="btnradio">Normal</label>"""
                    innerhtml += """<input id="toggle-off""" + str(
                        i) + """" class="toggle toggle-right" name="toggle""" + str(
                        i) + """" value="Perturbation" type="radio">"""
                    innerhtml += """<label for="toggle-off""" + str(
                        i) + """" class="btnradio">Perturbation</label></td></tr>"""
                i += 1
            innerhtml += """</table><br><button type="submit" class="btn" name="submit" id="analyze" value="analyze"><i class="fas fa-search"></i> Analyze</button>"""
            innerhtml += """<input type="hidden" id="hiddenfilename" name="hiddenfilename" value=\"""" + filename + """\"></form>"""
            html = """
            <script type="text/javascript">
            """
            html += """document.getElementById('datasetconfig').innerHTML = '""" + innerhtml + """';"""
            html += """
            </script>
            """
        return render_template('/dataset.html') + html


@server.route('/result', methods=['POST', 'GET'])
def result():
    if request.method == 'POST':
        global filename
        filename = request.form.get('hiddenfilename')
        df_headerdataset = gd.read_headerdataset(filename)
        sample = []
        genehead = ''
        i = 0
        for header in list(df_headerdataset):
            if i != 0:
                sample.append(request.form.get('toggle' + str(i)))
            else:
                genehead = header
            i += 1
        gd.assign_geneheadersampletype(genehead, sample)
        update = True
        if request.form.get('settinguser') == 'no':
            update = False
        html = ""
        datatable = gd.create_datatable(filename, update)
        datatablevalues = gd.create_datatablevalues(filename, update)
        bar = gd.create_plot(filename, update)
        library = gd.create_plotlibrary(filename, update)
        histogramcutoff = gd.create_histogramcutoff(filename, update)
        pca = gd.create_PCA(filename, update)
        cluster = gd.create_clustergrammer(filename, update)
        gd.rintegrate_test(filename, update)
        pathwaydb = gd.get_pathwaydb()
        barpathway = gd.gene_enrichment(filename, update)
        pathwayinfo = gd.pathwayinfo(filename, update)
        wordcloud = gd.create_wordcloud(filename)
        tabletargetup = gd.create_datatable_target('Up')
        tabletargetdown = gd.create_datatable_target('Down')
        valuestarget = gd.create_tabletarget(filename, update)
        dicttarget = gd.dict_target(filename)

        return render_template('/result.html', \
                               filename=filename, \
                               table=datatable, \
                               values=datatablevalues, \
                               plot=bar, \
                               plotlibrary=library, \
                               plothistogramcutoff=histogramcutoff, \
                               pcaplot=pca, \
                               clustergrammerplot=cluster, \
                               pathwaydb=pathwaydb, \
                               barpathwayplot=barpathway, \
                               pathwayinfo=pathwayinfo, \
                               wordcloud=wordcloud, \
                               tabletargetup=tabletargetup, \
                               tabletargetdown=tabletargetdown, \
                               valuestarget=valuestarget, \
                               dicttarget=dicttarget) + html


# @server.route('/fulldataset', methods=['POST', 'GET'])
# def fulldataset():
#     if request.method == 'POST':
#         dashp.continuedash(filename)
#     return redirect("http://localhost:5000/dash")
#
# dashp.initdash(server)

@server.route('/steps/stepone', methods=['POST', 'GET'])
def stepone():
    if request.method == 'POST':
        return render_template('/steps/stepone.html')

@server.route('/steps/steptwo', methods=['POST', 'GET'])
def steptwo():
    if request.method == 'POST':
        file = request.files['upload']
        filename = secure_filename(file.filename)
        file.save(os.path.join(server.config['UPLOAD_FOLDER'], filename))

        df_headerdataset = gd.read_headerdataset(filename)
        i = 0
        genehead = ''
        for header in list(df_headerdataset):
            if i == 0:
                genehead = header
                break
        gd.assign_geneheader(genehead)

        datatable = gd.create_datatable(filename, True)
        datatablevalues = gd.create_datatablevalues(filename, True)
        print(df_headerdataset)
        innerhtml = """"""
        innerhtml += """<table style="text-align: center;"><tr><th>Sample Name</th><th>Sample Type</th></tr>"""
        i = 0
        for header in list(df_headerdataset):
            if i != 0:
                innerhtml += """<tr><td>""" + header + """</td><td>"""
                innerhtml += """<input id="toggle-on""" + str(
                    i) + """" class="toggle toggle-left" name="toggle""" + str(
                    i) + """" value="Normal" type="radio" checked>"""
                innerhtml += """<label for="toggle-on""" + str(i) + """" class="btnradio">Normal</label>"""
                innerhtml += """<input id="toggle-off""" + str(
                    i) + """" class="toggle toggle-right" name="toggle""" + str(
                    i) + """" value="Perturbation" type="radio">"""
                innerhtml += """<label for="toggle-off""" + str(
                    i) + """" class="btnradio">Perturbation</label></td></tr>"""
            i += 1
        innerhtml += """</table>"""
        innerhtml += """<input type="hidden" id="hiddenfilename" name="hiddenfilename" value=\"""" + filename + """\">"""
        html = """
                    <script type="text/javascript">
                    """
        html += """document.getElementById('datasetconfig').innerHTML = '""" + innerhtml + """';"""
        html += """
                    </script>
        """

        return render_template('/steps/steptwo.html', table=datatable, values=datatablevalues) + html

@server.route('/steps/stepthree', methods=['POST', 'GET'])
def stepthree():
    if request.method == 'POST':

        filename = request.form.get('hiddenfilename')
        df_headerdataset = gd.read_headerdataset(filename)
        sample = []
        genehead = ''
        i = 0
        for header in list(df_headerdataset):
            if i != 0:
                sample.append(request.form.get('toggle' + str(i)))
            else:
                genehead = header
            i += 1
        gd.assign_geneheadersampletype(genehead, sample)

        library = gd.create_plotlibrary(filename, True)
        histogramcutoff = gd.create_histogramcutoff(filename, True)

        return render_template('/steps/stepthree.html', \
                               plotlibrary=library, \
                               plothistogramcutoff=histogramcutoff)

@server.route('/compoundanalysis', methods=['POST', 'GET'])
def compoundanalysis():
    if request.method == 'POST':
        target = request.form.get('hiddentarget')
        targetname = request.form.get('hiddentargetname')
        print(targetname)
        filename = request.form.get('hiddenfilenametarget')
        response = cm.getBioActs(filename, target)
        targetdetail = cm.gettargetdetail(target)

        if response == 'ok':
            response = ''
            bioact = cm.getgraphbioact(filename, target)
            allgraph = cm.getgraphpreprocess(filename, target)
            graphro5 = cm.RO5(filename, target)
            scatterro5 = cm.getscatterro5(filename, target)
            tablecompound = cm.create_datatablecompound(filename, target)
            datacompound = cm.create_datacompound(filename, target)
            cm.padel(filename, target)
            cm.modeling(filename, target)
            confustionplot = cm.getconfusion(filename, target)
            statvalue = cm.getstatvalue(filename, target)
            roccurve = cm.getplotroccurve(filename, target)
            elimination = cm.getelimination(filename, target)
            importance = cm.getimportance(filename, target)

            return render_template('/compoundanalysis.html', \
                                   filename=filename, \
                                   target=target, \
                                   targetname=targetname, \
                                   response=response, \
                                   targetdetail=targetdetail, \
                                   bioact=bioact, \
                                   allgraph=allgraph, \
                                   graphro5=graphro5, \
                                   scatterro5=scatterro5, \
                                   tablecompound=tablecompound, \
                                   datacompound=datacompound, \
                                   confustionplot=confustionplot, \
                                   statvalue=statvalue, \
                                   roccurve=roccurve, \
                                   elimination=elimination, \
                                   importance=importance)
        elif response == "Target doesn't not exist in Chembl." or response == "Target has no associated bioactivity or binding type assay.":
            return render_template('/compoundanalysis.html', \
                                   target=target, \
                                   targetname=targetname, \
                                   response=response, \
                                   targetdetail=targetdetail)

        bioact = cm.getgraphbioact(filename, target)
        return render_template('/compoundanalysis.html', \
                                   target=target, \
                                   targetname=targetname, \
                                   response=response, \
                                   targetdetail=targetdetail, \
                                   bioact=bioact)

@server.route('/aboutus', methods=['POST', 'GET'])
def aboutus():
    if request.method == 'POST':
        return render_template('/aboutus.html')

@server.route('/compoundinfo')
def compoundinfo():
    compound = request.args.get('compound', 0, type=str)
    moleculestructure = cm.getcompoundinfo(compound)
    return jsonify(result=moleculestructure)

@server.route('/compoundpredict', methods=['POST', 'GET'])
def compoundpredict():
    if request.method == 'POST':
        targetunidi = request.form.get('hiddentarget')
        filename = request.form.get('hiddenfilenametarget')
        file = request.files['upload']
        compoundfile = secure_filename(file.filename)
        path = os.path.join(server.config['UPLOAD_FOLDER'], filename + "_" + targetunidi)
        if not os.path.exists(os.path.join(path, 'compoundstest')):
            os.makedirs(os.path.join(path, 'compoundstest'))
        pathfilecompounds = os.path.join(path, 'compoundstest')
        file.save(os.path.join(pathfilecompounds, compoundfile))
        prediction = cm.predict(filename, targetunidi, compoundfile)
        return render_template('/compoundpredict.html', result=prediction)


if __name__ == '__main__':
    server.run(host='localhost', debug=True)
    # server.run(host='0.0.0.0')