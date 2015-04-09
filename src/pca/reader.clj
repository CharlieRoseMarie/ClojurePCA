(ns pca.reader)
(require '[clojure.data.csv :as csv]
         '[clojure.java.io :as io])

(def plastics "C:\\Users\\Tyler\\Projects\\sors\\clojure\\pca\\resources\\plastic_4_10.txt")

(defn read-file [filename & skip]
  (with-open [in-file (io/reader filename)]
    (doall
     (csv/read-csv in-file :separator \tab))))

(read-file plastics)
