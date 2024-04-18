const express = require('express')
const mod_fs = require('fs')
const app = express()
const puerto = 5000
const path = require('path')

const { exec } = require('child_process');

// const downloadPath = path.join(__dirname, 'uploads', '/')
// const fileInput = document.getElementById("download")

app.use(express.static(__dirname))

app.use(express.urlencoded({extended: true}))
app.use(express.json())

app.get("/", (req, res)=>{
    res.sendFile("index.html")
})

app.post("/datosform", (req, res)=>{

    let datos_html = false 

    try {
        datos_html = mod_fs.readFileSync("./datos.html", {encoding: 'utf-8', flag: 'r'})
        datos_html = datos_html.replace("%nombre%", req.body.nombre)
        datos_html = datos_html.replace("%apellido%", req.body.apellido)
    } catch (error) {
        res.status(200).send("Archivo no encontrado. " + error)
        
    }
    return res.send(datos_html)
})

app.get("/formconGet", (req, res)=>{
    let nombre = req.query.nombreg
    let apellido = req.query.apellidog

    if (nombre.length <= 0 || apellido.length <= 0) {
        return res.status(200).json({mensaje: "Debe de capturar nombre y apellido"})   
    }
    return res.status(200).json({mensaje: "Datos del formulario con GET", datos: 
    {nombre: nombre, apellido: apellido}})
})

app.listen(puerto, ()=>{
    console.log("servidor iniciado");
})


app.get("/descargar/:id", (req, res)=>{
    res.download(__dirname+'/uploads/'+req.params.id, req.params.id,(err)=>{
        if(err) {
            console.log(err)
            res.status(404).send("no such file or directory")
            
        }else{
            console.log('LISTO')
        }
    })

})


app.get('/python_script/:function/:nameFig',async (req, res)=>{
    let file_down = false 

    try {
        file_down = mod_fs.readFileSync("./filedown.html", {encoding: 'utf-8', flag: 'r'})
        // datos_html = datos_html.replace("%nombre%", req.body.nombre)
        // datos_html = datos_html.replace("%apellido%", req.body.apellido)
    } catch (error) {
        res.status(200).send("Archivo no encontrado. " + error)
        
    }


    try{
        const fun = req.params.function;
        const namefig = req.params.nameFig;
        // const namefig = `${Date.now()}-${req.params.nameFig}`;
        // console.log(`function: ${fun}, namefig: ${namefig}`);

        exec(`python child-processes/script.py ${fun} ${namefig}`, (error, stdout, stderr) => {
        // exec(`python src/uploads/script_python3.py ${fun} ${namefig}`, (error, stdout, stderr) => {
            if (error) {
            console.error(`error: ${error.message}`);
            res.send(error.message)
            return;
            }
        
            if (stderr) {
            console.error(`stderr: ${stderr}`);
            res.send(stderr)

            return;
            }
        
            console.log(`stdout:\n${stdout}`);
            res.status(200).send(file_down)
            // res.status(200).send(`${output_script}`)

        })
    } catch (error) {
        res.status(500).json(error);
    }
          
} );
